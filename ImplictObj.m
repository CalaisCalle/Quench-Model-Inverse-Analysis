function [errout] = ImplictObj(h_input_vec)
%IMPLICTOBJ Summary of this function goes here
%   Detailed explanation goes here
load('params.mat');

%set the h values of the parameters 
params.hs = h_input_vec; %h_bot

thermo_data = importdata(params.fname, '\t', 1);
data_time = thermo_data.data(:,1);
data_temp = thermo_data.data(:,2:end);

%% Create Prediction Data
[model_time,model_temp] = ImplicitModel(params);
% Interpolate Temperatures
[n_row, n_col] = size(model_temp);

% interpolate times
min_data_time = min(data_time);
min_model_time = min(model_time);
min_time = max([min_data_time, min_model_time]);
max_time = min([max(model_time), max(data_time)]);
inter_time = min_time: max_time / 150: max_time;

error_sum = 0; % Variable to hold error
for i = 1:n_col
   %interpolate data and model
   inter_model_T = interp1(model_time, model_temp(:,i), inter_time);
   inter_data_T = interp1(data_time, data_temp(:,i), inter_time);
   
   %Calculate Error on column
   adj = 1E-10; % Adjustment value to avoid dividing by zero
   location_error = abs(inter_data_T - inter_model_T) ./ (inter_data_T + adj);
   error_sum = error_sum + trapz(inter_time, location_error);
end

t_max = inter_time(end);
errout = error_sum / t_max;

end

