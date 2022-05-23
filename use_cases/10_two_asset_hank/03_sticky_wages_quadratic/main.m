%------------------------------------------------------------------------%
% 
% This code computes a two-asset HANK model with sticky wages.
% 
% Code written by Andreas Schaab.
% Current version: May 2022. First version: May 2022.
% 
% If you find this code helpful in your own work, please cite:
%  1. Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%     with Adaptive Sparse Grids. Working Paper.
%  2. Dávila, E. and A. Schaab. Optimal Monetary Policy with Heterogeneous
%     Agents: A Timeless Ramsey Approach. Working Paper.
% Thanks!
% 
%------------------------------------------------------------------------%

clear
close all
clc
warning off

diary ./output/output.log
diary on

addpath(genpath('../../../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS
param = define_parameters('shock_type', 'demand', 'T', 120, 'N', 100, 'shock_theta', log(2), ...
                          'implicit_g', 1, 'xi', 0, 'rho', 0.08, 'delta', 0.03, ...
                          ... % 'Delta_KF', 0.05, 'maxit_KF', 100000, ...
                          'tau_lab', 0);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'k'});
G_dense.dx = G_dense.da * G_dense.dk;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'k'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);

% Initialize BCs:
for j = 1:param.discrete_types
    BC{1}.left.type = '0'; BC{1}.right.type = '0';
    BC{2}.left.type = '0'; BC{2}.right.type = '0';
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end

% For convenience:
G.aa = [G.a, G.a];
G.kk = [G.k, G.k];
G.zz = [repmat(param.zz(1), [G.J, 1]), repmat(param.zz(2), [G.J, 1])];


%% 0-INFLATION STATIONARY EQUILIBRIUM
fprintf('\n\n:::::   0-INFLATION STATIONARY EQUILIBRIUM   :::::: \n\n');

piw = 0;

% Get better guess for value function:
r0 = 0.80 * param.rho; K0 = 6; N0 = 0.8; X0 = [r0, K0, N0];
[~, G, G_dense, ~] = stationary(X0, piw, G, G_dense, param);

% Solve for steady state prices:
options = optimset('Display', 'iter', 'UseParallel', false, 'TolX', 1e-15);
X = fsolve(@(x) stationary(x, piw, G, G_dense, param), X0, options);

% Solve with correct prices:
[~, G, G_dense, ss] = stationary(X, piw, G, G_dense, param);

fprintf('Stationary Equilibrium:  r = %.4f   K = %.4f   N = %.4f \n', ss.r, ss.K, ss.N);
fprintf('Markets:  goods=%.1d   bonds=%.1d   labor=%.1d   capital=%.1d   savings=%.1d \n\n', ...
    ss.excess_goods, ss.excess_bonds, ss.excess_labor, ss.excess_capital, ss.excess_saving);


%% COMPUTE TRANSITION DYNAMICS
fprintf('\n\n:::::::::::   TRANSITION DYNAMICS   ::::::::::: \n\n');

% Shock:
switch param.shock_type
    case 'TFP'
        z0 = ss.Z * ones(param.N, 1);
    case 'demand'
        z0 = param.rho * ones(param.N, 1);
    case 'cost-push'
        z0 = param.epsilon * ones(param.N, 1);
end
dz = param.shock_level * ones(param.N, 1);
for n = 1:param.N-1
    dz(n+1, :) = exp(-param.shock_theta * param.t(n+1)) .* param.shock_level;
end

z = z0 + dz;

% Initialize policy and guesses:
t0 = zeros(param.N, 1);
Y0 = ss.Y * ones(param.N, 1); 
I0 = ss.I * ones(param.N, 1);
M0 = ss.M * ones(param.N, 1);
x0 = [Y0; I0; M0];

diff0 = transition(x0, t0, z, ss, G, G_dense, param, 'markets');

% Solve for price paths:
f = @(x, y) transition(x, t0, z, ss, y{1}, y{2}, param, 'markets'); y0{1} = G; y0{2} = G_dense;
x = fsolve_newton(f, reshape(x0, [numel(x0), 1]), diff0, y0, 0, 5, 2);
sim = transition(x, t0, z, ss, G, G_dense, param, 'all');


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

% DSS: 
% for n = 1:adapt_iter
%     
%     figure('visible', 'off');
%     l1=scatter3(G_adapt{n}.a, G_adapt{n}.k, V_adapt{n}(:, 1));
%     xlabel('Liquid: $a$', 'Interpreter', 'Latex');
%         xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
%         set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
%     ylabel('Illiquid: $k$', 'Interpreter', 'Latex');
%         ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
%         set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
%     zlabel('$V^U(a, k)$', 'Interpreter', 'Latex');
%         zlh = get(gca, 'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
%         set(zlh, 'Rotation', 0, 'Position', zlp+[0, 11, -2], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');        
%     set(gcf, 'renderer', 'Painters');
%     exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);
% 
% end


% Transition:
figure('visible', 'on');
subplot(3, 3, 1); 
plot(param.t, 100 * (sim.Y - ss.Y)/ss.Y); ylabel('% dev'); title('$Y_t$', 'Interpreter', 'Latex');
subplot(3, 3, 2); 
plot(param.t, 100 * (sim.K - ss.K)/ss.K); title('$K_t$', 'Interpreter', 'Latex');
subplot(3, 3, 3); 
plot(param.t, 100 * (sim.N - ss.N)/ss.N); title('$N_t$', 'Interpreter', 'Latex');
subplot(3, 3, 4); 
plot(param.t, 100 * (sim.C - ss.C)/ss.C); ylabel('% dev'); title('$C_t$', 'Interpreter', 'Latex');
subplot(3, 3, 5); 
plot(param.t, 100 * (sim.I - ss.I)/ss.I); title('$I_t$', 'Interpreter', 'Latex');
subplot(3, 3, 6); 
plot(param.t, 100 * (sim.w - ss.w)/ss.w); title('$w_t$', 'Interpreter', 'Latex');
subplot(3, 3, 7); 
plot(param.t, 100 * (sim.r - ss.r)); ylabel('%'); xlabel('Quarters'); title('$r_t$', 'Interpreter', 'Latex');
subplot(3, 3, 8); 
plot(param.t, 100 * (sim.rk - ss.rk)); xlabel('Quarters'); title('$r_t^k$', 'Interpreter', 'Latex');
subplot(3, 3, 9); 
plot(param.t, z); xlabel('Quarters'); title([param.shock_type, ' shock (lvl)'], 'Interpreter', 'Latex');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/transition_dynamics.eps');


diary off


