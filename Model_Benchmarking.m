%% Benchmarking
% Using the objective function to benchmark the code against data

%% set parameters
params.rho_As = [7840 0.44];
params.k_As = [13.1947 0.0126919];
params.Cp_As = [490 0.0733333];

params.T0 = 1000;
params.T_inf = 293;
params.eps = 1;
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
ns = [10,20,30,40,50];
% array to save errors
errors = zeros(length(ns),1);
times = zeros(length(ns),1);

% loop through n
for i = 1:length(ns)
    params.n = ns(i);
    save('params.mat', 'params');
    tic
    errors(i) = ImplictObj(hs);
    times(i) = toc;
end

figure(1)
subplot(1,2,1)
plot(ns, errors*100,'-x')

figure(1)
subplot(1,2,2)
plot(ns, times, '-x')