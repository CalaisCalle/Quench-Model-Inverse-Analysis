function [errout] = ImplictObj(h_input_vec)
%IMPLICTOBJ Summary of this function goes here
%   Detailed explanation goes here

load('params.mat', params);
params.hs = h_input_vec;

[t,predict_mat] = ImplicitSolver(params);

fname = 'thermocouple_data.txt';

thermo_data = importdata(fname, '\t', 1);
thermo_t = thermo_data.data(:,1);
thermo_temp = thermo_data.data(:,2:10);

T = predict_mat;
T_intp = interp1(t,T,thermo_t);

err = abs(thermo_temp - T_intp)./thermo_temp;

for j = 1: 9
    % Dunno if the loop is necessary tbh
    interr = trapz(thermo_t,err(:,j));
    errout = errout + interr;
end
t_max = thermo_t(end);
errout = errout / t_max;
end

