params.rho_As = [4430 0];
params.k_As = [1.117 0.0174];
params.Cp_As = [546.31 0.219];

params.fname = 'thermocouple_data.txt';
params.T0 = 960 + 273;
params.T_inf = 293;
params.n = 20;
params.eps = 0.279;
params.CFL = 0.1;
params.time = 0;
params.t_max = 18000;
params.L = 0.25;
params.xs = [0.5, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5, 0.05, 0.05]*params.L; % X Points
params.ys = [0.5, 0.05, 0.05, 0.05, 0.5, 0.95, 0.95, 0.95, 0.5]*params.L; % Y Points
save('params.mat', 'params');
% hs = [50 100 20];

h0(1) = 30;
h0(2) = 30;
h0(3) = 30;

h_LB(1) = 10;
h_LB(2) = 10;
h_LB(3) = 10;

h_UB(1) = 150;
h_UB(2) = 150;
h_UB(3) = 150;
tic
[hs, err] = fminsearchbnd('ImplictObj',h0,h_LB,h_UB);
fprintf('Time taken: %fs\n', toc);
fprintf('h_top: %f\nh_sides: %f\nh_bot: %f\nerror: %f%', hs, err*100);
