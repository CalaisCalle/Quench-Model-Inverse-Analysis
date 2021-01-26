%% Benchmarking
% Using the objective function to benchmark the code against data
% Used to verify the model and choose parameters.

%% set parameters
params.rho_As = [7840 0.44];
params.k_As = [13.1947 0.0126919];
params.Cp_As = [490 0.0733333];

params.T0 = 1000;
params.T_inf = 293;
params.eps = 1;
params.n = 20;
params.CFL = 0.1;
params.time = 0;
params.L = 0.03;
% Thermistor points
params.xs = [0.015 0.027 0.027]; % X Points
params.ys = [0.015 0.015 0.027]; % Y Points

params.fname = 'Benchmark_data.txt';

% Set HTCs
hs = [50, 100, 20];
% Range of n to benchmark over
ns = 10:10:100;
CFLs = 1:-0.1:0.1;
% array to save errors
errors = zeros(length(ns),1);
times = zeros(length(ns),1);
%% Fixed CFL
%%  loop through n
for i = 1:length(ns)
    params.n = ns(i);
    save('params.mat', 'params');
    %% get execution time and error
    tic
    errors(i) = ImplictObj(hs);
    times(i) = toc;
end

%% Plot n against time and error 
% error
figure(1)
subplot(1,2,1)
plot(ns, errors*100,'-x')
xlabel('n')
ylabel('Error %')

figure(1)
subplot(1,2,2)
plot(ns, times, '-x')
xlabel('n')
ylabel('Execution Time [s]')
saveas(gcf, 'n_error_time.pdf')

%% Fixed n, change CFL
% reset errors and set n
errors = zeros(length(CFLs), 1);
times = errors;
params.n = 20;

for i = 1:length(CFLs)
    params.CFL = CFLs(i);
    save('params.mat', 'params');
    tic
    errors(i) = ImplictObj(hs); % reuse same hs
    times(i) = toc;
end

%% Plots for CFL study
% error
figure(2)
subplot(1,2,1)
plot(CFLs, errors*100, '-x')
xlabel('CFL')
ylabel('Error %')

% Execution time
figure(2)
subplot(1,2,2)
plot(CFLs, times, '-x')
xlabel('CFL')
ylabel('Execution Time [s]')
% Save plot 2
saveas(gcf, 'cfl_error_time.pdf')
