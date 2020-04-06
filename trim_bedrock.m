% trim out flow trajectory that crash into bedrock
clc
clear all
close all

load ('flowpath_MC.mat')

%just need to grab the 72m sample, correspond to 8th element in 1st
%dimension of h_age_save and x_age_save

h_age_72 = h_age_save(8,:,:);
x_age_72 = x_age_save(8,:,:);

h_age_72 = reshape (h_age_72, [861 1000]);
x_age_72 = reshape (x_age_72, [861 1000]);

%find the NaNs
%in reality the bedrock extend to 72100 m upstream of the glacier. Ice
%trajectory that crash into the bedrock, at 72000m upstream if we
%interpolate it will be NaN

for i = 1:1000
    xx = [72100-100];
    x_age = x_age_72(:,i);
    [x_age, index] = unique(x_age); 
    h_age = h_age_72(:,i);
    interp(i) = interp1 (x_age,h_age(index),xx);
end

index = isnan(interp);

x_trim = x_age_save(:,:,~index);
h_trim = h_age_save(:,:,~index);

%Calculate weight probability for the flow
flow_rand = flow_rand(~index);
P_n_I = exp(-0.5.*(flow_rand).^2);

save ('flowpath_trim.mat', 'x_trim', 'h_trim','P_n_I','flow_rand','age')

%optional: try to plot the surviving trajectories
x_plot = reshape (x_trim (8,:,:), [861 989]);
h_plot = reshape (h_trim (8,:,:), [861 989]);

figure (1)
plot (x_plot, Elev_age+h_plot)
hold on
plot (FL_x, (FL_bed))


