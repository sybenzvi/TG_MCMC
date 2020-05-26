# Markov Chain Monte Carlo Algorithm for Taylor Glacier 14C data

Forked version from Michael Dyonisius.

This readme file presents basic descriptions for the Markov Chain Monte Carlo (MCMC) algorithm written in MATLAB. The MCMC method is used to constrain two parameters in cosmogenic 14C production model by muons. For simplicity, the two model parameters that are optimized with the MCMC method are &quot;_fneg_&quot; and &quot;_ffast_&quot; – corresponding scaling factors for negative muon capture and fast muon reaction that are constant with respect to depths. The MCMC method aims to optimize these 2 parameters given the observations – which in this case represent 14 total 14C measurements from 7 unique depth levels (6.85m, 15m, 19.5m, 40m, 51m, 61.5m, 72m) obtained from ice cores drilled in Taylor Glacier, Antarctica.

## File descriptions

1. **External data**<br>
  a. flowpath\_MC.mat – contains the pool of 1000 flowpaths (see section 1)<br>
  b. flowpath\_trim.mat – contains the flowpaths that didn&#39;t crash into bedrock (see section 2)<br>
  c. P\_neutron.mat – production rate from neutron as function of depth (see section 3)<br>
  d. P\_muon.mat – production rate from muons as function of depth (see section 3)<br>
  e. all\_data.mat – contains measurement data, sample depths and uncertainties<br>
2. **Functions**<br>
  a. trim\_trajectory.m – this function removes part of the flow trajectory that goes out of bounds (goes beyond 72km), see section 2<br>
  b. calc14C\_v3.m – the main integration function for 14C concentration in the ice parcel (see section 3)<br>
  c. interp1qr_m – function for faster 1D interpolation<br>
3. **Run scripts**<br>
  a. trim\_bedrock.m – trim flowpath\_MC.mat into flowpath\_trim.mat (see section 2)<br>
  b. run\_MCMC.m – run the MCMC algorithm (see section 4)<br>

## 1. Ice flow model

Buizert et al. (2012) developed a 2D ice flow model for Taylor Glacier. Using estimates of ablation rates along the glacier and bedrock profile from Kavanaugh et al. (2009a) the model generates an ice parcel back-trajectory for each given sample depth. For the MCMC algorithm I perturbed the ablation rates within their uncertainties and generated 1000 flow trajectories for 8 unique depth levels (6.85m, 10m, 15m, 19.5m, 40m, 51m, 61.5m, 72m). I added the 10m sample depth because we are only missing 14CH4 data on that depth level, and 14CH4 account for less than 1% of total 14C signal. At the moment, the 10m data are not included in the MCMC fit algorithm, but later they can be if necessary. The 8x1000 flow trajectories are saved in &quot; **flowpath\_MC.mat**&quot; file.

Important variables in the &quot; **flowpath\_MC.mat**&quot; file include

- _age_ (1x861) – this is the age of the ice parcel trajectory, going from 0 at the drill location to 6000 years before present.
- _h\_age\_save_ (8x861x1000) – this is the depth of the ice parcel trajectory (y-location in the glacier). The first dimension (8) represents the 8 unique depth levels. The 2nd dimension (861) corresponds to the _age_ and the 2nd dimension of _x\_age\_save_ (x-location in the glacier). The 3rd dimension (1000) represents the pool of 1000 possible ice parcel trajectories given the perturbation in ablation rates.
-  _x\_age\_save_ (8x861x1000) – this is the x-location of the ice parcel trajectory. The first dimension (8) represents the 8 unique depth levels. The 2nd dimension (861) corresponds to the _age_ and the 2nd dimension of _h\_age\_save_ (y-location in the glacier). The 3rd dimension (1000) represents the pool of 1000 possible ice parcel trajectories given the perturbation in ablation rates.
- _flow\_rand_ (1x1000) – this is the random number (η) that corresponds to the likelihood of each 1000 possible ice parcel trajectories. The probability for a given _flow\_rand_ is
<br>**P(η|I) = exp(-0.5\*((η- η****best****)/σ)^2)**<br> where ηbest = mean of η and σ = stdev of η. because the _flow\_rand_ variable is generated with MATLAB&#39;s randn function, ηbest = mean = 0 and σ = stdev = 1, thus P(η|I) simplifies into
<br>**P(η|I) = exp(-0.5\*(flow\_rand)^2)**

## 2. Trimming the trajectories.

Under extreme ablation rate scenarios, the ice parcel trajectories crash into bedrock and become unphysical afterwards. I use the **&quot;trim\_bedrock.m&quot;** run script to simply remove the trajectories that crash into bedrock and resave the remaining trajectories in **&quot;flowpath\_trim.mat&quot;** file.

Additionally, the bedrock data only extend to 72km upstream of the glacier terminus. Given different ablation rates scenario, some ice parcel trajectories go &quot;out of bounds&quot; beyond 72km faster than others and the flow trajectories also become unphysical afterwards since there&#39;s no bedrock.

To account for the limitation in the extent of the bedrock data, we assume that the amount of 14C at the beginning of the trajectory (corresponding to 72km in the x-direction) is at steady state (production equals decay). The age of our 72m deep sample is approximately 70 kyr BP. Survey by Morse et al. (1998) showed that the depth of ~70 kyr BP ice at northern flank of Taylor Dome is ~575m, thus we assumed that the depth of the long term transport for the 72m sample is 575m. At 72km, the 72m sample under the &quot;best&quot; ablation scenario is at 699m deep, thus we use 699m as the baseline depth &quot;z\_baseline.&quot; For other ice parcel trajectories, the steady state depth of long-term transport (z\_deep) is calculated following

**z\_deep = 575 – (z\_baseline – zx=72km)**

where zx=72km is the depth of the corresponding sample at 72km x-direction.

We want to trim the unphysical parts out so that the integration of 14C content in the model can be consistent. MATLAB doesn&#39;t like having a matrix with different number of elements, thus the function **&quot;trim\_trajectory.m&quot;** takes an individual ice parcel trajectory [_h\_age (1x861), x\_age(1x861), age(1x861),z\_baseline(1x1)_] and transform them into cell {} format where it is possible to save individual vectors with different elements in them.

## 3. 14C production model

The baseline 14C production rate by neutrons, in unit of 14C atoms g ice-1 yr-1 as function of depth (m) is obtained from LSD model by Lifton et al. (2014) and saved in **&quot;P\_neutron.mat&quot;** file. The baseline 14C production rate by muons, in unit of 14C atoms g ice-1 yr-1 as function of depth (m) is obtained from LSD model by Balco et al. (2008) and saved in **&quot;P\_muon.mat&quot;** file. Both baseline 14C production rate files are scaled to the altitude and latitude of Taylor Glacier drill site.

From the baseline production rates, the amount of 14C in the sample for a given flow trajectory is calculated by the differential equation
**dC/dt = fneg\*Pneg(z) + ffast\*Pfast(z) – λC**
fneg         = constant scaling factor for negative muon capture (unitless)
Pneg         = 14C production from negative muon capture (14C atoms g ice-1 yr-1)
fneg         = constant scaling factor for fast muon reaction (unitless)
Pneg         = 14C production from fast muon reaction (14C atoms g ice-1 yr-1)
z        = sample depth (m) – note that z is also a function of time z(t)
λ         = radiocarbon decay constant (1.216e-4 yr-1)
C         = 14C concentration in the ice parcel (14C atoms g ice-1)

As an initial condition (C0) before the integration, we assumed that at the depth of long-term transport (z\_deep), the 14C concentration in the ice parcel is at equilibrium
dC/dt (at z\_deep) = 0 = fneg\*P\_neg\*z\_deep + ffast\* P\_fast\*z\_deep - C0\*λ
currently we are ignoring production from neutron because we&#39;re only fitting to data below 6.85m.

In the MATLAB code, the expected total 14C in the ice parcel at the end of flow trajectory (_Cexp_) is calculated with the function **&quot;calc14C\_v3.m&quot;.** The **&quot;calc14C\_v3.m&quot;** function takes _P\_n (production rate from neutron), P\_neg (production rate from negative muon), P\_muf (production rate from fast muon), z\_P (depth corresponding to the production rate), h\_age (depth vector of the ice trajectory), age (age vector of the age trajectory), C0 (initial condition), and f (3x1 matrix that linearly scale the production rates, [fn fneg ffast])._

## 4. Markov Chain Monte Carlo algorithm**

The MCMC is run with the **&quot;run\_MCMC.m&quot;** script. The general recipe for the MCMC is

1. Start with guess values of fneg and ffast
2. Pick 100 random trajectory scenarios among the 1000 possibilities
3. Calculate the expected 14C (_Cexp_) in the samples for the 100 trajectory scenarios
4. Calculate the probability of observed 14C concentration {Cobs} for given {Cexp}<br>
**P({Cobs}|{Cexp}) = prod (1/sqrt(2\*pi)\*(1/{Cobs\_error})\*exp(-0.5\*(({Cobs} – {Cexp})/{Cobs\_error)^2))**<br>
_{Cexp}                = vector containing all expected 14C concentrations_<br>
_{Cobs}                = vector containing measured 14C_<br>
_{Cobs\_error}        = vector containing analytical uncertainty in measured 14C_<br>
calculate P({Cobs}|{Cexp}) for all 100 flowlines, each corresponding to η (_flow\_rand_)<br>
5. Numerically calculate the integrated probability<br>
**Integrate dη P({Cobs}|{Cexp}) P(η|I)**<br>
In MATLAB, I simply used the trapz function (trapezoidal numerical integration) to numerically integrate P({Cobs}|{Cexp}) P(η|I) vs. η
6. Now generate new values of fneg and ffast<br>
_fneg\_new = f\_neg\_old + randn(1) \* spacing parameter_<br>
_ffast\_new = f\_fast\_old + randn(1) \* spacing parameter_<br>
The spacing parameter determines how far the fneg and ffast move around. Currently it is set as 0.01
7. Repeat step 2-5 with _fneg\_new_ and _ffast\_new_
8. Calculate the ratio R = Pnew/Pold
9. If R>= uniform random number between 0 to 1<br>
then move to a new location (accept the new fneg and ffast)<br>
else<br>
stay with the old fneg and ffast<br>
