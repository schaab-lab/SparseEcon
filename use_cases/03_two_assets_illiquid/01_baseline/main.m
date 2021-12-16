%------------------------------------------------------------------------%
% 
% This code solves a two-asset portfolio choice problem with one liquid
% (bond) and one illiquid (capital) asset using adaptive sparse grids. 
% The model is a simplified version of the partial equilibrium household 
% problem in Schaab (2020).
% 
% Code written by Andreas Schaab.
% Current version: December 2021. First version: December 2019.
% 
% If you find this code helpful in your own work, please cite:
%  1. Schaab, A. Micro and Macro Uncertainty. Working Paper.
%  2. Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%     with Adaptive Sparse Grids. Working Paper.
% Thanks!
% 
% This code replicates the DSS of HANK6 (in Schaab, 2020) for the first 
% iteration before grid adaptation:
% [r = 0.0309, K=5.1977, L=0.7404, Q=1.0000]   [G=1.9e-08, B=1.3e-09, L=-6.3e-10, C=-6.5e-09])
% After correcting param.L definition (and associated KF lines):
% [r = 0.0308, K=5.2309, L=0.7430, Q=1.0000]   [G=1.8e-08, B=-3.2e-10, L=7.7e-11, C=2.2e-09]
% 
%------------------------------------------------------------------------%

clear
close all
clc

diary ./output/output.log
diary on

addpath(genpath('../../../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS

param = define_parameters('add_tol', 1e-4, 'keep_tol', 5e-6, 'max_adapt_iter', 3);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a','k'});
G_dense.dx = G_dense.da*G_dense.dk;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a','k'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);

% Initialize BCs:
for j = 1:param.discrete_types
    BC{1}.left.type = '0'; BC{1}.right.type = '0';
    BC{2}.left.type = '0'; BC{2}.right.type = '0';
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.002; K0 = 7; L0 = 0.8; Q0 = 1; tau0 = 0.3;
    if exist('ss', 'var'), r0 = ss.r; K0 = ss.K; L0 = ss.L; Q0 = ss.Q; tau0 = ss.tau; end
    X0 = [r0, K0, L0, Q0, tau0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    % options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    % X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf(['Stationary Equilibrium: [r = %.4f, K=%.4f, L=%.4f, Q=%.4f]   ', ...
             'Markets=[G=%.1d, B=%.1d, L=%.1d, C=%.1d])\n\n'],...
              ss.r, ss.K, ss.L, ss.Q, ss.excess_goods, ss.excess_bonds, ss.excess_labor, ss.excess_capital);
    
    V_adapt{adapt_iter} = ss.V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
    
    % Update BCs:
    for j = 1:param.discrete_types
        BC{1}.left.type = '0'; BC{1}.right.type = '0';
        BC{2}.left.type = '0'; BC{2}.right.type = '0';
        G = gen_FD(G, BC, num2str(j));
    end
end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off');
    l1 = scatter3(G_adapt{n}.a, G_adapt{n}.k, V_adapt{n}(:, 1));
    xlabel('Liquid: $a$', 'Interpreter', 'Latex');
        xlh = get(gca,'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Illiquid: $k$', 'Interpreter', 'Latex');
        ylh = get(gca,'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    zlabel('$V^U(a,k)$', 'Interpreter', 'Latex');
        zlh = get(gca,'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
        set(zlh, 'Rotation', 0, 'Position', zlp+[0, 11, -2], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');        
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off


