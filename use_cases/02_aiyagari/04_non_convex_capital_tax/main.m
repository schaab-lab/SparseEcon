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

% profile('-memory','on');
% ps = parallel.Settings;
% ps.Pool.AutoCreate = false;

diary ./output/output.log
diary on

addpath(genpath('../../../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS

param = define_parameters();


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
    if stats.n_change == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
    
end


%% FINE GRID + COMPARISONs
G_fine = setup_grid(0, param.surplus_fine, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
G_fine.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_fine);

% GFine.V0 = sparseProject(ss.V, GFine.grid, G);
[~, ~, ~, ss_fine] = stationary([ss.K, ss.tau], G_fine, G_dense, param);

% G.BHFine = get_projection_matrix(GFine.grid, GFine.lvl, G);

% fprintf('L_1   norm |V - Vfine| = %.2d\n', norm(G.BHFine * ss.V - ssFine.V, 1));
% fprintf('L_2   norm |V - Vfine| = %.2d\n', norm(G.BHFine * ss.V - ssFine.V, 2));
% fprintf('L_inf norm |V - Vfine| = %.2d\n', norm(G.BHFine * ss.V - ssFine.V, Inf));

% G.normInf = norm(G.BHFine * ss.V - ssFine.V, Inf);
% clear GFine.BHDense G.BHFine;

% Sampling points:
P = 30000; points = rand(P, 2);
G_fine.V_test = zeros(P,1);

for k = 1:P
    G_fine.V_test(k) = sparse_project(ss_fine.V, points(k, :), G_fine);
end

% Uniform dense grids:
for j = 1:param.fine_level % j = max corresponds to "fine" grid
    
    G_tpg{j} = setup_grid(0, j-1+param.surplus, param.min, param.max, ...
        'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
    G_tpg{j}.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_tpg{j});
    
    [~, ~, ~, ss_tpg{j}] = stationary([ss.K, ss.tau], G_tpg{j}, G_dense, param);
    
    % GTPG{j}.BHFine = get_projection_matrix(GFine.grid, GFine.lvl, GTPG{j});
    
    % GTPG{j}.normInf = norm(GTPG{j}.BHFine * ssTPG{j}.V - ssFine.V, Inf);
    % norm(sparseProject(ssTPG{j}.V, GFine.grid, GTPG{j}) - ssFine.V, Inf);
    % clear GTPG{j}.BHFine GTPG{j}.BHDense;
    
    G_tpg{j}.V_test = zeros(P, 1);
    for k = 1:P
        G_tpg{j}.V_test(k) = sparse_project(ss_tpg{j}.V, points(k, :), G_tpg{j});
    end
    G_tpg{j}.norm_inf = norm(G_tpg{j}.V_test - G_fine.V_test, Inf);
    
    fprintf('L_inf norm |V_tpg%.i - V_fine| = %.2d  with # grid points: %.i\n', ...
        j, G_tpg{j}.norm_inf, G_tpg{j}.J);

end


% Regular sparse grids:
for j = 1:param.fine_level+3
    
    G_sg{j} = setup_grid(j-1, param.surplus, param.min, param.max, ...
        'NamedDims', {1, 2}, 'Names', {'k', 'z'}, 'DxxDims', 2);
    G_sg{j}.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G_sg{j});
    
    % GSG{j}.V0 = sparse_project(ss.V, GSG{j}.grid, G);
    [~, ~, ~, ss_sg{j}] = stationary([ss.K, ss.tau], G_sg{j}, G_dense, param);
    
    % GSG{j}.BHFine = get_projection_matrix(GFine.grid, GFine.lvl, GSG{j});    
    % GSG{j}.normInf = norm(GSG{j}.BHFine * ssSG{j}.V - ssFine.V, Inf);
    % clear GSG{j}.BHFine GSG{j}.BHDense;
    
    G_sg{j}.V_test = zeros(P,1);
    for k = 1:P
        G_sg{j}.V_test(k) = sparse_project(ss_sg{j}.V, points(k,:), G_sg{j});
    end
    G_sg{j}.norm_inf = norm(G_sg{j}.V_test - G_fine.V_test, Inf);

    fprintf('L_inf norm |V_sg%.i - Vfine| = %.2d  with # grid points: %.i\n', ...
        j, G_sg{j}.norm_inf, G_sg{j}.J);

end


% Adaptive sparse grids:
for j = 1:adapt_iter
    
    % G_adapt{j}.BHFine  = get_projection_matrix(GFine.grid, GFine.lvl, G_adapt{j});
    % clear G_adapt{j}.BHFine;
    
    G_adapt{j}.V_test = zeros(P,1);
    for k = 1:P
        G_adapt{j}.V_test(k) = sparse_project(G_adapt{j}.V0, points(k, :), G_adapt{j});
    end
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

figure;
subplot(1, 2, 1); hold on;
l1 = plot(1:numel(norm_asg), norm_asg);
l2 = plot(1:numel(norm_sg), norm_sg);
l3 = plot(1:numel(norm_tpg), norm_tpg);
hold off; title('$L_\infty$ norm','Interpreter','Latex');
xlabel('Grid levels / adaptations');
legend([l1, l2, l3],{'Adaptive Sparse Grid', 'Regular Sparse Grid', 'Uniform Dense Grid'}, 'box', 'off');

subplot(1, 2, 2); hold on;
l1 = plot(1:numel(grid_points_asg), grid_points_asg);
l2 = plot(1:numel(grid_points_sg), grid_points_sg);
l3 = plot(1:numel(grid_points_tpg), grid_points_tpg);
hold off; title('Grid points', 'Interpreter', 'Latex');
xlabel('Grid levels / adaptations');
    
set(gcf,'renderer','Painters');
exportgraphics(gcf, './output/error1.eps');


figure; hold on;
l1 = plot(grid_points_asg, norm_asg); 
l2 = plot(grid_points_sg,  norm_sg);
l3 = plot(grid_points_tpg, norm_tpg);
hold off;
xlabel('Grid points'); 
title('$L_\infty$ norm', 'Interpreter', 'Latex');
set(gca, 'XScale', 'log', 'YScale', 'log');
    
legend([l1,l2,l3],{'Adaptive Sparse Grid', 'Regular Sparse Grid', 'Uniform Dense Grid'}, 'box', 'off');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/error2.eps');


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off');
    subplot(1, 2, 1);
    l1 = scatter3(G_adapt{n}.k, G_adapt{n}.z, V_adapt{n});
    xlabel('Capital: $k$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh,'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    title('$V(a,z)$', 'Interpreter', 'Latex');
        
    subplot(1, 2, 2);
    l1=scatter3(G_adapt{n}.k, G_adapt{n}.z, c_adapt{n});
    xlabel('Capital: $k$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$','Interpreter','Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment','middle', 'HorizontalAlignment', 'right');
    title('$c(a,z)$', 'Interpreter', 'Latex');
    
    set(gcf, 'renderer', 'Painters', 'Position', [10 10 1200 400]);
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end


diary off

% profile off
% profile report


