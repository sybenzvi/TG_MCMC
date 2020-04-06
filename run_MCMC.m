%run MCMC algorithm to optimize fneg and ffast

%clear all previously loaded data
clc
clear all
close all

%load necessary external data
load ('P_muon.mat')
load ('P_neutron.mat')
load ('flowpath_trim.mat')
load('all_data.mat')
C14_error = C14_error./2; %the error in .mat data are 2sigma so divide by 2
f_best = [0,0.205,0.25]; %this should be the best answer for the fractions, from the "bruteforce" method
%1st column is for neutron, 2nd column is for negative muon, 3rd column is
%for fast muon

P_n = interp1(Z_neutron,Lal_P_n,z_P); %the production rate by neutron (P_neutron) is at different depth resolution
%than the muons, thus interpolate so they are consistent.

z_baseline = -699.3415; %this is 72m deep sample, at 72km mark, on baseline trajectory

% Constants:
lambda = 1.216e-4;  %14C decay constant, yr^-1
ro = 0.92;          % ice density, g/cm^3

repeat_data = [1 1 3 3 4 4 5 5 6 6 7 7 8 8]; %repeat depth vector

%convert to column vector for faster interpolation
z_P = z_P.';
P_muf = P_muf.';
P_n = P_n.';
P_neg = P_neg.';

%% 1st Markov Chain run
omega = 100; %number of flowpaths per random set of f's

f_new = [0,0.3,0.3]; %starting guess

  %pick uniform index of random flow
    in = round(rand (omega,1).*(989-1)+1); %uniform number between 1 to 989 
    
    %trim_trajectory that goes out of bounds
    for i = 1:omega
        %first put the x and y of the flowlines into temp variables
        temp_x = x_trim(:,:,in(i));
        temp_h = h_trim(:,:,in(i));

        temp_x = reshape (temp_x,[8 1 861]);
        temp_h = reshape (temp_h,[8 1 861]);
        
        %then calculate depth of long term transport
        z_deep = 575 -( temp_h(:,861) - z_baseline);
    
        
        for i2 = 1:8 %i2 = number of unique sample depths
            h_age_temp = temp_h(i2,:,:);
            x_age_temp = temp_x(i2,:,:);
            [x_flow{i2},y_flow{i2},age_flow{i2},z_deep(i2,:)] = trim_trajectory(h_age_temp,x_age_temp,age,z_baseline); %variables are saved as cells so it can take different amount of indexes
        end
        
        %now the cell {x_flow}, {y_flow}, {age_flow} should have 8 elements
        %in them, each element should contain 1x790ish elements
        %corresponding to the trimmed flow

        parfor i3 = 1:8 %i3 = also number of unique sample depth but for C14 integration
        %calculate initial condition  
            P0_neg = interp1qr(z_P,P_neg,z_deep(i3));
            P0_muf = interp1qr(z_P,P_muf,z_deep(i3));
            P0 = f_new(2).*P0_neg+f_new(3).*P0_muf;
            C0 = P0./lambda;
            %calculate expected 14C given fneg, ffast, and flow that we
            %just trim above
            Cexp(i3,:) = calc14C_v3(P_n,P_neg,P_muf,z_P,y_flow{i3},age_flow{i3},C0,f_new);
        end
        
        Cexp_repeat = Cexp(repeat_data); %repeat Cexp so that the vector corresponds to the observations (including replicates)
        
        %calculate P(Cobs|Cexp) for all 14 observations vs. expected  
        P(i,:) = 1/sqrt(2*pi).*(1./C14_error).*exp(-0.5.*((C14_sample- Cexp_repeat)./C14_error).^2);        
    end

%Calculate P({Cobs}|{Cexp}) by multiplying all P(Cobs|Cexp) in the same rows     
P_sum = prod (P,2);
n_weight = P_n_I(in).'; %P(n|I) was already calculated in trim_bedrock.m so just pick P(n|I) correspond to index in

%P({Cobs}|{Cexp}) P(?|I)
prob_product = P_sum.*n_weight;

%remove NaNs so numerical integration works
nan_index = isnan (prob_product);
prob_product(nan_index) = [];

%do the same with the eta (n) parameter
n = flow_rand(in).';
n (nan_index) = [];

%combine P({Cobs}|{Cexp}) P(?|I) and n into single matrix
A = sortrows ([n prob_product]);
%Integrate d? P({Cobs}|{Cexp}) P(?|I)
new_P = trapz (A(:,1),A(:,2));

Psave (1,:) = new_P;

%figure (1) %uncomment this to plot P({Cobs}|{Cexp}) P(?|I) vs. n 
%plot (A(:,1),A(:,2));

%% Establish logical counter
move = true;
f_save = f_new;
old_P = new_P;
f_old = f_new;

%% Run markov chain
chain_length = 1000; %markov chain length

for i4 = 2:chain_length

    if move == true
        f2_new = f_old(2) + randn(1).*0.01;
        f3_new = f_old(3) + randn(1).*0.01;
        f_new = [0 f2_new f3_new];
    end
    
    in = round(rand (omega,1).*(989-1)+1);
    
    for i = 1:omega
        temp_x = x_trim(:,:,in(i));
        temp_h = h_trim(:,:,in(i));

        temp_x = reshape (temp_x,[8 1 861]);
        temp_h = reshape (temp_h,[8 1 861]);
        z_deep = 575 -( temp_h(:,861) - z_baseline);
    
        %trim_trajectory that goes out of bounds
        for i2 = 1:8 %i2 = number of unique sample depths
            h_age_temp = temp_h(i2,:,:);
            x_age_temp = temp_x(i2,:,:);
            [x_flow{i2},y_flow{i2},age_flow{i2},z_deep(i2,:)] = trim_trajectory(h_age_temp,x_age_temp,age,z_baseline); %variables are saved as cells so it can take different amount of indexes
        end

        parfor i3 = 1:8 %i3 = also number of unique sample depth but for C14 integration
        %initial condition  
            P0_neg = interp1qr(z_P,P_neg,z_deep(i3));
            P0_muf = interp1qr(z_P,P_muf,z_deep(i3));
            P0 = f_new(2).*P0_neg+f_new(3).*P0_muf;
            C0 = P0./lambda;
            Cexp(i3,:) = calc14C_v3(P_n,P_neg,P_muf,z_P,y_flow{i3},age_flow{i3},C0,f_new);
        end
    
        Cexp_repeat = Cexp(repeat_data);
        P(i,:) = 1/sqrt(2*pi).*(1./C14_error).*exp(-0.5.*((C14_sample- Cexp_repeat)./C14_error).^2);
        
    end

    P_sum = prod (P,2);
    n_weight = P_n_I(in).';
    prob_product = P_sum.*n_weight;
    
%remove NaNs so numerical integration works
nan_index = isnan (prob_product);
prob_product(nan_index) = [];

%do the same with the eta (n) parameter
n = flow_rand(in).';
n (nan_index) = [];

%combine P({Cobs}|{Cexp}) P(?|I) and n into single matrix
A = sortrows ([n prob_product]);
%Integrate d? P({Cobs}|{Cexp}) P(?|I)
new_P = trapz (A(:,1),A(:,2));
    
%rejection criteria
    u = rand;
    r = new_P./old_P;
    
    if u <= r
        move = true;
        old_P = new_P;
        fprintf('move true ')
    else
        move = false;
        fprintf('move false ')
    end
    
    movetracker (i4) = move;
    rsave (i4) = r;
    usave (i4) = u;
    Psave (i4,:) = new_P;
    f_save(i4,:) = f_new;
    fprintf('run #%d\n', i4)
end

f_all = f_save (:,(2:3));
%save ('temp_results.mat')
return
