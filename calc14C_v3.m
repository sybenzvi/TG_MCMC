function [C_fin] = calc14C_v3(P_n,P_neg,P_muf,z_P,h_age,age,C0,f)
% function integrates over flowpath with euler's method

age_temp = cell2mat(age); %temporary age, need to do it because ageflow is cell
tstep = 0.5;
age_interp = [age_temp(1):tstep:age_temp(end)]; %generate age vector to integrate
    
%age 0 = initial parcel at depth because integration needs to be forward
%and positive NOTE: this is the opposite of Christo's age
    [age_temp, index] = unique(cell2mat(age)); %need to trim out non-unique ages
    y_temp = cell2mat(h_age);
    y_temp = y_temp(index);
    parcel_depth = fliplr(interp1(age_temp,y_temp,age_interp).*-1); %generate parcel depth history corresponding to age interp
    C(1) = C0;

 for i = 1:(length(age_interp)-1)
        
        Pn_interp = f(1).*interp1qr(z_P,P_n,parcel_depth(i));
        Pneg_interp = f(2).*interp1qr(z_P,P_neg,parcel_depth(i));
        Pmuf_interp = f(3).*interp1qr(z_P,P_muf,parcel_depth(i));
                
    dCdt = Pn_interp+Pneg_interp+Pmuf_interp-C(i)*1.216e-4;
    C(i+1) = C(i)+tstep.*dCdt;
 end
 
C_fin = C(end);