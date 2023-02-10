%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a model in the spirit 
%   of Cagetti and De Nardi (2006) and Buera and Shin (2013)
%   using adaptive sparse grids.
% The contiunous time formulation of the model is taken
%   from https://benjaminmoll.com/codes/ 
% 
% Code written by Sergi Barcons and Andreas Schaab.
%   Current version: December 2022. First version: December 2022. 
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

addpath(genpath('../../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS

param = define_parameters('add_tol',1e-4,'keep_tol',1e-5);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);
G_dense.dx = G_dense.da*G_dense.dz;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID

blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.045; if exist('ss', 'var'), r0 = ss.r; end; X0 = r0; %J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 2);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: (r = %.4f, w = %.4f),  markets(K = %.2f,  S = %.2d,  Y - δK - C = %.2d) \n\n', ...
        ss.r, ss.w, round(ss.K_s,2), ss.S, ss.excess_supply);
    
    V_adapt{adapt_iter} = ss.V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);

end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off');
    l1 = scatter3(G_adapt{n}.a, G_adapt{n}.z, V_adapt{n}(:, 1));
    %xticks([param.amin, 5, 10, 15, param.amax]); xticklabels({num2str(param.amin), '', '', '', num2str(param.amax)});
    %xticklabels({});
    xlabel('Wealth: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Ability: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh,'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    zlabel('$V(a,z)$', 'Interpreter', 'Latex');
        zlh = get(gca, 'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
        set(zlh, 'Rotation', 0, 'Position', zlp+[0, 0.08, 0], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off




