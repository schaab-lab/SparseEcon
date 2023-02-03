%------------------------------------------------------------------------%
% 
% PREAMBLE
% 
%------------------------------------------------------------------------%

clear
close all
clc

diary ./output/output.log
diary on

addpath(genpath('../../../lib/'))
addpath(genpath('/Users/andreasschaab/Dropbox/numericalMethods/lib/export_fig'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS
param = define_parameters('T', 15, 'N', 150, 'shock_type', 'demand', 'shock_theta', log(2));


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(param.l_dense, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'a'});
G_dense.dx = G_dense.da;

% Sparse grid:
G = setup_grid(param.l, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'a'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);

% For convenience:
G.aa = [G.a, G.a];
G.zz = [repmat(param.zz(1), [G.J, 1]), repmat(param.zz(2), [G.J, 1])];


%% STATIONARY EQUILIBRIUM
fprintf('\n\n:::::   STATIONARY EQUILIBRIUM   :::::: \n\n');

piw = 0;

% Get better guess for value function:
r0 = 0.80 * param.rho; N0 = 1; X0 = [r0, N0];
[~, G, G_dense, ~] = stationary(X0, G, G_dense, param);

% Solve for steady state prices:
options = optimset('Display', 'iter', 'UseParallel', false, 'TolX', 1e-12);
X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);

% Solve with correct prices:
[~, G, G_dense, ss] = stationary(X, G, G_dense, param);

fprintf('Stationary Equilibrium: [r = %.4f, N = %.2f],  markets[S=%.2d,  Y-C=%.2d, B=%.2d, Union=%.2d] \n\n', ...
    ss.r, ss.N, ss.S, ss.excess_supply, ss.excess_bonds, ss.excess_union);

% FD matrix:
Da_dense{1} = (G_dense.DS_interior.('D1F'){1} + G_dense.(['DS_', '1']).('D1F'){1}) .* (G.BH_dense * ss.s(:, 1) > 0) + ...
              (G_dense.DS_interior.('D1B'){1} + G_dense.(['DS_', '1']).('D1B'){1}) .* (G.BH_dense * ss.s(:, 1) < 0) + ...
              (G_dense.DS_interior.('D1C'){1} + G_dense.(['DS_', '1']).('D1C'){1}) .* (G.BH_dense * ss.s(:, 1) == 0);
Da_dense{2} = (G_dense.DS_interior.('D1F'){1} + G_dense.(['DS_', '2']).('D1F'){1}) .* (G.BH_dense * ss.s(:, 2) > 0) + ...
              (G_dense.DS_interior.('D1B'){1} + G_dense.(['DS_', '2']).('D1B'){1}) .* (G.BH_dense * ss.s(:, 2) < 0) + ...
              (G_dense.DS_interior.('D1C'){1} + G_dense.(['DS_', '2']).('D1C'){1}) .* (G.BH_dense * ss.s(:, 2) == 0);

Da_dense{1}(1, 1) = abs(Da_dense{1}(2, 2));
Da_dense{1}(1, 2) = 0;

if ss.s(1, 2) <= 0, error('Employed household is not saving at borrowing constraint.\n'); end

if ss.s(end, 1) >= 0
    Da_dense{1}(end, end) = abs(Da_dense{1}(end-1, end-1));
    Da_dense{1}(end, end-1) = -abs(Da_dense{1}(end-1, end-1));
end
if ss.s(end, 2) >= 0
    Da_dense{2}(end, end) = abs(Da_dense{2}(end-1, end-1));
    Da_dense{2}(end, end-1) = -abs(Da_dense{2}(end-1, end-1));
end

ss.Da_dense = Da_dense;
Da = blkdiag(Da_dense{1}, Da_dense{2});
ss.Da = Da;


%% TRANSITION DYNAMICS: GLOBAL
fprintf('\n\n::::::     TRANSITION DYNAMICS: GLOBAL     :::::: \n\n');

% Shock:
switch param.shock_type
    case 'TFP'
        z0 = ss.Z * ones(param.N, 1);
    case 'demand'
        z0 = param.rho * ones(param.N, 1);
    case 'cost-push'
        z0 = param.epsilon * ones(param.N, 1);
    case 'monetary'
        z0 = zeros(param.N, 1);
end

dz = param.shock_level * ones(param.N, 1);
for n = 1:param.N-1
    dz(n+1, :) = exp(-param.shock_theta * param.t(n+1)) .* param.shock_level;
end

z = z0 + dz;

% Initialize guesses:
Y0 = ss.Y * ones(param.N, 1); M0 = ss.M * ones(param.N, 1);
x0 = [Y0; M0];

% Transition dynamics: global
diff0 = transition(x0, z, ss, G, G_dense, param, 'markets');
f = @(x, y) transition(x, z, ss, y{1}, y{2}, param, 'markets'); y0{1} = G; y0{2} = G_dense;
x = fsolve_newton(f, reshape(x0, [numel(x0), 1]), diff0, y0, 0, 5, 2);
sim{1} = transition(x, z, ss, G, G_dense, param, 'all');

% run_irfs(sim, ss, param);


%% TRANSITION DYNAMICS: DIRECT
fprintf('\n\n::::::     TRANSITION DYNAMICS: DIRECT     :::::: \n\n');

query = {'H_x', 'H_z'};
f = @(x, z) transition(x, z, ss, G, G_dense, param, 'markets');
H = get_jacobians_market_clearing(x0, z0, f, param, query, 1);

dx = - inv(H.H_x) * H.H_z * dz;

x = x0 + dx;
sim{2} = transition(x, z, ss, G, G_dense, param, 'all');

run_irfs(sim, ss, param);


%% TRANSITION DYNAMICS: FAKE NEWS
fprintf('\n\n::::::     TRANSITION DYNAMICS: FAKE NEWS     :::::: \n\n');

H_direct = H; clear H;

run_time = tic;
% [H, p] = fake_news_testing(x0, z0, ss, G, G_dense, param);
[H, p] = fake_news(x0, z0, ss, G, G_dense, param);
run_time = toc(run_time); fprintf('Fake-news algorithm run-time: %.2f seconds\n', run_time);

fprintf('Max difference in H_z : %.2d\n', max(max(abs(H.H_z - H_direct.H_z))));
fprintf('Max difference in H_x : %.2d\n', max(max(abs(H.H_x - H_direct.H_x))));

dx = - inv(H.H_x) * H.H_z * dz;

x = x0 + dx;
sim{3} = transition(x, z, ss, G, G_dense, param, 'all');

run_irfs(sim, ss, param);


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

% Optimal Policy IRF:
fprintf('\nPlotting Figures...\n');

diary off






