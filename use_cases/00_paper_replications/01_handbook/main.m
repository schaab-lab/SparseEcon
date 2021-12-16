%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a model similar to
% Aiyagari (1994) using adaptive sparse grids. Labor supply is inelastic 
% and earnings risk follows a diffusion process. We introduce a non-linear
% capital income tax tau(k) = 0 for k<k* and tau(k) = kappa for k>k*.
% 
% Code written by Johannes Brumm, Christopher Krause, Andreas Schaab and
% Simon Scheidegger.
% Current version: December 2021. First version: October 2021.
% 
% If you find this code helpful in your own work, please cite:
%  1. Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%     with Adaptive Sparse Grids. Working Paper.
%  2. Brumm, J., C. Krause, A. Schaab, and S. Scheidegger. Sparse Grids for 
%     Dynamic Economic Models. 2021. Available at SSRN 3979412.
% Thanks!
% 
% The code for our Handbook Chapter is at:
% https://github.com/SparseGridsForDynamicEcon/SparseGrids_in_econ_handbook
% 
%------------------------------------------------------------------------%

clear
close all
clc
warning off

% profile('-memory', 'on');
% ps = parallel.Settings;
% ps.Pool.AutoCreate = false;

diary ./output/output.log
diary on

addpath(genpath('../../../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS

% param = define_parameters('surplus', [4, 1], 'fine_level', 8, 'l_dense', [7, 2], ...
%                           'kmin', 0, 'kmax', 50, ...
%                           'max_adapt_iter', 8, 'kappa', 0.04, 'kstar', 27);
param = define_parameters('surplus', [5, 2], 'fine_level', 7, 'l_dense', [7, 2], ...
                          'kmin', 0, 'kmax', 50, ...
                          'rho', 0.02, 'theta_z', 0.10, 'sig_z', 0.04, ...
                          'max_adapt_iter', 7, 'kappa', 0.00, 'kstar', 40);

%{ 
    Observations:
    - lDense does not seem to make a big / any difference
    - If stationary doesn't run on refined GSG, lower Delta (~0.5 works)
    -
%}

%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
G_dense.dx = G_dense.dk * G_dense.dz;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1); 
c_adapt = cell(param.max_adapt_iter, 1); 
s_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    K0 = 10; tau0 = 0; if exist('ss', 'var'), K0 = ss.K; tau0 = ss.tau; end; X0 = [K0, tau0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: (r = %.4f, K = %.2f),  markets(S = %.2d,  Y-C-I = %.2d, Kgap = %.2d) \n\n', ...
        ss.r, ss.K, ss.S, ss.excess_supply, ss.excess_capital);
    
    V_adapt{adapt_iter} = ss.V; c_adapt{adapt_iter} = ss.c; s_adapt{adapt_iter} = ss.s; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    if adapt_iter == param.max_adapt_iter, break; end
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.nChange == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
    
end


%% FINE GRID + COMPARISONs

% Sampling points:
rng(1); P = 30000; points = rand(P, 2);


% Fine grid:
% G_fine = setup_grid(0, param.surplus_fine+1, param.min, param.max, ...
%     'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
G_fine = setup_grid(0, param.surplus_fine, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
G_fine.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_fine);

% G_fine.V0 = sparse_project(ss.V, G_fine.grid, G);
[~, ~, ~, ss_fine] = stationary([ss.K, ss.tau], G_fine, G_dense, param);

V_test = zeros(P, 1);
for k = 1:P
    V_test(k) = sparse_project(ss_fine.V, points(k, :), G_fine);
end
G_fine.V_test = V_test;

% Uniform dense grids:
for j = 1:param.fine_level % j = max corresponds to "fine" grid
    
    G_tpg{j} = setup_grid(0, j-1+param.surplus, param.min, param.max, ...
        'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
    G_tpg{j}.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_tpg{j});
    
    [~, ~, ~, ss_tpg{j}] = stationary([ss.K, ss.tau], G_tpg{j}, G_dense, param);
    
    V_test = zeros(P, 1);
    for k = 1:P
        V_test(k) = sparse_project(ss_tpg{j}.V, points(k, :), G_tpg{j});
    end
    G_tpg{j}.V_test = V_test;
    G_tpg{j}.norm_inf = norm(G_tpg{j}.V_test - G_fine.V_test, Inf);
    
    fprintf('L_inf norm |V_tpg%.i - V_fine| = %.2d  with # grid points: %.i\n', ...
        j, G_tpg{j}.norm_inf, G_tpg{j}.J);

end


% Regular sparse grids:
% CAUTION: if this breaks, use param.Delta=0.5; param.maxit=2000;
for j = 1:param.fine_level %+2
    
    G_sg{j} = setup_grid(j-1, param.surplus, param.min, param.max, ...
        'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
    G_sg{j}.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_sg{j});
    
    [~, ~, ~, ss_sg{j}] = stationary([ss.K, ss.tau], G_sg{j}, G_dense, param);
    
    V_test = zeros(P, 1);
    for k = 1:P
        V_test(k) = sparse_project(ss_sg{j}.V, points(k, :), G_sg{j});
    end
    G_sg{j}.V_test = V_test;
    G_sg{j}.norm_inf = norm(G_sg{j}.V_test - G_fine.V_test, Inf);

    fprintf('L_inf norm |V_sg%.i - V_fine| = %.2d  with # grid points: %.i\n', ...
        j, G_sg{j}.norm_inf, G_sg{j}.J);

end


% Adaptive sparse grids:
for j = 1:adapt_iter
    
    [~, ~, ~, ss_asg{j}] = stationary([ss.K, ss.tau], G_adapt{j}, G_dense, param);

    V_test = zeros(P, 1);
    for k = 1:P
        V_test(k) = sparse_project(ss_asg{j}.V, points(k, :), G_adapt{j});
    end
    G_adapt{j}.V_test = V_test;
    G_adapt{j}.norm_inf = norm(G_adapt{j}.V_test - G_fine.V_test, Inf);
    
end


% Plot:
norm_asg = []; norm_sg = []; norm_tpg = []; grid_points_asg = []; grid_points_sg = []; grid_points_tpg = [];
for j = 1:max([numel(G_adapt), numel(G_sg), numel(G_tpg)]) %param.fine_level-1
    if j > numel(G_adapt)
        continue;
    else
        norm_asg(j) = G_adapt{j}.norm_inf;
        grid_points_asg(j) = G_adapt{j}.J;
    end
    
    if j > numel(G_sg)
        continue;
    else
        norm_sg(j)  = G_sg{j}.norm_inf;
        grid_points_sg(j)  = G_sg{j}.J;
    end
    
    if j > numel(G_tpg)
        continue;
    else
        norm_tpg(j) = G_tpg{j}.norm_inf;
        grid_points_tpg(j) = G_tpg{j}.J;
    end    
end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off');
    subplot(1, 2, 1);
    l1 = scatter3(G_adapt{n}.k, G_adapt{n}.z, V_adapt{n});
    xlabel('Capital: $k$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    title('$V(a, z)$', 'Interpreter', 'Latex');
        
    subplot(1, 2, 2);
    l1=scatter3(G_adapt{n}.k, G_adapt{n}.z, c_adapt{n});
    xlabel('Capital: $k$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    title('$c(a, z)$', 'Interpreter', 'Latex');
    
    set(gcf, 'renderer', 'Painters', 'Position', [10 10 1200 400]);
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

figure('visible', 'off');
subplot(1, 2, 1); hold on;
l1 = plot(1:numel(norm_asg), norm_asg);
l2 = plot(1:numel(norm_sg), norm_sg);
l3 = plot(1:numel(norm_tpg), norm_tpg);
hold off; title('$L_\infty$ norm', 'Interpreter', 'Latex');
xlabel('Grid levels / adaptations');
legend([l1, l2, l3], {'Adaptive Sparse Grid', 'Regular Sparse Grid', 'Uniform Dense Grid'}, 'box', 'off');

subplot(1, 2, 2); hold on;
l1 = plot(1:numel(grid_points_asg), grid_points_asg);
l2 = plot(1:numel(grid_points_sg), grid_points_sg);
l3 = plot(1:numel(grid_points_tpg), grid_points_tpg);
hold off; title('Grid points', 'Interpreter', 'Latex');
xlabel('Grid levels / adaptations');
    
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/error1.eps');


figure('visible', 'off'); hold on;
l1 = plot(grid_points_asg, norm_asg); 
l2 = plot(grid_points_sg,  norm_sg);
l3 = plot(grid_points_tpg, norm_tpg);
hold off;
xlabel('Grid points'); 
title('$L_\infty$ norm', 'Interpreter', 'Latex');
set(gca, 'XScale', 'log', 'YScale', 'log');
    
legend([l1, l2, l3], {'Adaptive Sparse Grid', 'Regular Sparse Grid', 'Uniform Dense Grid'}, 'box', 'off');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/error2.eps');

diary off

% profile off
% profile report



%% PLOTS FOR HANDBOOK
n = adapt_iter;

figure('visible', 'off');
scatter3(G_adapt{n}.k, G_adapt{n}.z, V_adapt{n}); 
ylim([0.3, 1.5]);
xlabel('Capital: $k$', 'Interpreter', 'Latex', 'FontSize', 24);
    xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
    set(xlh, 'Rotation', 14, 'Position', [1, 0.8, 1].*xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
ylabel('Earnings: $z$', 'Interpreter', 'Latex', 'FontSize', 24);
    ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', -25, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
title('$V(k, z)$', 'Interpreter', 'Latex', 'FontSize', 30);

set(gcf, 'renderer', 'Painters', 'Position', [10 10 600 400]);
if param.kappa > 0
    exportgraphics(gcf, './output/handbook_aiyagari_nonconvex_full_adaptation_V.eps');
else
    exportgraphics(gcf, './output/handbook_aiyagari_vanilla_full_adaptation_V.eps');
end

figure('visible', 'off');
scatter3(G_adapt{n}.k, G_adapt{n}.z, c_adapt{n});
ylim([0.3, 1.5]);
xlabel('Capital: $k$', 'Interpreter', 'Latex', 'FontSize', 24);
    xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
    set(xlh, 'Rotation', 14, 'Position', [1, 0.8, 1].*xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
ylabel('Earnings: $z$', 'Interpreter', 'Latex', 'FontSize', 24);
    ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', -25, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
title('$c(k, z)$', 'Interpreter', 'Latex', 'FontSize', 30);

set(gcf, 'renderer', 'Painters', 'Position', [10 10 600 400]);
if param.kappa > 0
    exportgraphics(gcf, './output/handbook_aiyagari_nonconvex_full_adaptation_c.eps');
else
    exportgraphics(gcf, './output/handbook_aiyagari_vanilla_full_adaptation_c.eps');
end


figure('visible', 'off'); hold on;
l1 = plot(grid_points_asg, norm_asg, '-x', 'LineWidth', 1.5, 'MarkerSize', 8); 
l2 = plot(grid_points_sg,  norm_sg , '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
l3 = plot(grid_points_tpg, norm_tpg, '-*', 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;
xlabel('Grid points', 'FontSize', 20); 
title('$L_\infty$ norm', 'Interpreter', 'Latex', 'FontSize', 30);
ylim([max([norm_asg(n-1), norm_sg(n-1)]), max([norm_asg, norm_sg, norm_tpg])]);
set(gca, 'XScale', 'log', 'YScale', 'log');
    
legend([l1, l2, l3], {'ASG', 'SG', 'Full Grid'}, 'box', 'off', 'FontSize', 20);
set(gcf, 'renderer', 'Painters'); box on;
if param.kappa > 0
    exportgraphics(gcf, './output/handbook_aiyagari_nonconvex_error.eps');
else
    exportgraphics(gcf, './output/handbook_aiyagari_vanilla_error.eps');
end


