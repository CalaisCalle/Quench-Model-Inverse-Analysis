%% Set Material Parameters/Coefficients
params.rho_As = [4430 0]; % Density [Kg m^-3]
params.k_As = [1.117 0.0174]; % Thermal conductivity [W m^-1 k^-1] 
params.Cp_As = [546.31 0.219]; % Specific Heat Capacity [J K^-1 Kg^-1]
params.eps = 0.279; % Emissivity

params.fname = 'thermocouple_data.txt';
%% Set material dimensions and quantisation
params.L = 0.25; % length [m]
params.n = 20; % Grid points on one edge (total is n^2)

%% Initial Temperatures
params.T0 = 960 + 273; % Initial temperature [K]
params.T_inf = 293; % Far-Field Temperature [K]

%% Time and time step conditions
params.CFL = 0.1; % Courant-Friedrich-Lewy Constant
params.time = 0;
params.t_max = 18000;

%% Thermistor locations
params.xs = [0.5, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5, 0.05, 0.05]*params.L; % X Points
params.ys = [0.5, 0.05, 0.05, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5]*params.L; % Y Points
save('params.mat', 'params');
% hs = [50 100 20];

%% Initial values of heat transfer coefficients [W m^-2]
h0(1) = 30; % h_top
h0(2) = 30; % h_side
h0(3) = 30; % h_bot

%% Lower bound h values
h_LB(1) = 10; % h_top
h_LB(2) = 10; % h_side
h_LB(3) = 10; % h_bot

h_UB(1) = 150; % h_top
h_UB(2) = 150; % h_side
h_UB(3) = 150; % h_bot

tic % time it

%% Run optimisation save error to file
% There's a typo in the name of the object function - need to change
[hs, err] = fminsearchbnd('ImplictObj',h0,h_LB,h_UB);
result_file = fopen('Results.txt', 'w'); % open file
% print results to file
fprintf(result_file, 'Time taken: %fs\n', toc);
fprintf(result_file,'h_top: %f\nh_sides: %f\nh_bot: %f\n', hs);
fprintf(result_file,'Error: %.3f\n', err);
fclose(result_file);

%% Create graphs for comparison
% set h values to those found from optimisation
params.hs = hs;
% Run model to get final prediction
[model_time_f, model_data_f] = ImplicitModel(params);

% load data from file
thermo_data = importdata(params.fname, '\t', 1);
data_time = thermo_data.data(:,1);
data_temp = thermo_data.data(:,2:end);

% Create figure with subplots
[n_rows, n_cols] = size(model_data_f);
figure(1)
for i = 1:n_cols
    subplot(3,3,i) % Have to set 3x3 subplots manually: could improve
    % plot data as solid red line
    plot(model_time_f, model_data_f(:,i) - 273, 'kx', 'MarkerSize', 5);
    hold on
    % plot prediction as series of x
    plot(data_time, data_temp(:,i) -273, '-r'); % plot in celsius
    hold off
    xlabel('Time [s]')
    ylabel('Temperature [\circC]')
    legend('Data', 'Sim')
    title(['Position', num2str(i)]);

end
saveas(gcf, 'sim_model_comparison.pdf') % save the figure