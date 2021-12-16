%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a two-asset life-cycle 
% model similar to Kaplan-Violante (ECTA 2014) using adaptive sparse grids.
% 
% Code written by Andreas Schaab.
% Current version: December 2021. First version: July 2021.
% 
% If you find this code helpful in your own work, please cite:
%   Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%   with Adaptive Sparse Grids. Working Paper.
% Thanks!
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

param = define_parameters('add_tol', 5e-4, 'keep_tol', 5e-8, 'max_adapt_iter', 5, ...
                          'deathrate', 0, 'tau_lab', 0.15, 'UI', 0.15, 'gov_bond_supply', 0.20);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2, 3}, 'Names', {'a', 'k', 't'});
G_dense.dx = G_dense.da*G_dense.dk*G_dense.dt;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2, 3}, 'Names', {'a', 'k', 't'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    K0 = 5; L0 = 0.75; Q0 = 1; tau0 = 0.25;
    if exist('ss', 'var'), K0 = ss.K; L0 = ss.L; Q0 = ss.Q; tau0 = ss.tau; end
    X0 = [K0, L0, Q0, tau0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 2);
    options = optimset('Display', 'iter', 'UseParallel', true, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf(['Stationary Equilibrium: [r = %.4f, K=%.4f, L=%.4f, Q=%.4f]   ', ...
             'Markets=[G=%.1d, B=%.1d, L=%.1d, C=%.1d])\n\n'], ...
              ss.r, ss.K, ss.L, ss.Q, ss.excess_goods, ss.excess_bonds, ss.excess_labor, ss.excess_capital);
    
    V_adapt{adapt_iter} = ss.V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    if adapt_iter==param.max_adapt_iter, break; end
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change==0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);

end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off');
    l1=scatter3(G_adapt{n}.a, G_adapt{n}.t, V_adapt{n}(:, 1));
    xlabel('Liquid: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Age: $t$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    zlabel('$V^U(a, :, t)$', 'Interpreter', 'Latex');
        zlh = get(gca, 'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
        set(zlh, 'Rotation', 0, 'Position', zlp+[0, 11, -2], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');        
    set(gcf, 'renderer', 'Painters');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off


