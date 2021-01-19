function [errout] = ImplictObj(h_input_vec)
%IMPLICTOBJ Summary of this function goes here
%   Detailed explanation goes here
load('params.mat');

%% set the heat transfer coefficient for the solver from input 
params.hs = h_input_vec; 

%% Load the data and assign it 
thermo_data = importdata(params.fname, '\t', 1);
data_time = thermo_data.data(:,1);
data_temp = thermo_data.data(:,2:end);

%% Create Prediction Data
[model_time,model_temp] = ImplicitModel(params);

% Get number of rows/columns for the purpose of indexing data
% Doing it this way makes the code easier to adapt to different conditions
[n_row, n_col] = size(model_temp);

%% Created interpolated times
% Get the highest of the 2 minimum times and lowest of the 2 maximum times
min_time = max([min(data_time),  min(model_time)]);
max_time = min([max(model_time), max(data_time)]);
% Resulting times ensured to be within limits of both data and model times
inter_time = min_time: max_time / 150: max_time;

%% Calculate model error
error_sum = 0; % Variable to hold error
for i = 1:n_col
   %interpolate data and model
   inter_model_T = interp1(model_time, model_temp(:,i), inter_time);
   inter_data_T = interp1(data_time, data_temp(:,i), inter_time);
   
   %% Calculate Error on column of data
   adj = 1E-10; % Adjustment value to avoid dividing by zero
   location_error = abs(inter_data_T - inter_model_T) ./ (inter_data_T + adj);
   error_sum = error_sum + trapz(inter_time, location_error);
end

%% Return an average error
t_max = inter_time(end);
errout = error_sum / t_max;

end

