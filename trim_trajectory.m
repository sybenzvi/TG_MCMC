function [x_flow,y_flow,age_flow,z_init] = trim_trajectory(h_age,x_age,age,z_baseline)
%This fuction trims the ice parcel trajectories so they can be consistently
%integrated from 72km out to 14km out.
%z baseline is 72m deep sample, at 72km mark, on baseline trajectory

% reindex flow history into cell
h_age = reshape (h_age, [1 861]);
x_age = reshape (x_age, [1 861]);

%find index where the x flow goes out of bounds
[ix,iy] = find (x_age >= 72000);

%trim them out
x_age_trim = x_age;
h_age_trim = h_age;

x_age_trim(ix,iy) = NaN;
h_age_trim(ix,iy) = NaN;


x_age_temp = x_age_trim; %temp var 
h_age_temp = h_age_trim;
nan_index = isnan(x_age_temp);
    
%remove nans
x_age_temp = x_age_temp(~nan_index);
h_age_temp =h_age_temp(~nan_index);
age_temp = age(~nan_index);
    
%extend last data point to 72000 so that all ice parcel have consistent
%starting point
x_age_temp = [x_age_temp 72000];
    
[x_age_interp, index] = unique(x_age);
h_age_temp = [h_age_temp interp1(x_age_interp,h_age(index),72000)];
age_temp = [age_temp round(interp1(x_age_interp,age(index),72000))]; %do round otherwise integration doesnt work
    
    x_flow{1} = x_age_temp;
    age_flow{1} = age_temp;
    y_flow{1} = h_age_temp;
    z_init = 575 -( h_age_temp(end) - z_baseline);

end

