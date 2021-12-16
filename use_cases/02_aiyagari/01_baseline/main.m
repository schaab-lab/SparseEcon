%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a model similar to
% Aiyagari (1994) using adaptive sparse grids. Labor supply is endogenous 
% and earnings risk follows a two-state Markov chain.
% 
% Code written by Andreas Schaab and Allen Zhang.
% Current version: September 2021. First version: September 2019.
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

param = define_parameters();


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(param.l_dense, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'k'});
G_dense.dx = G_dense.dk;

% Sparse grid:
G = setup_grid(param.l, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'k'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    K0 = 30; L0 = 1; if exist('ss', 'var'), K0 = ss.K; L0 = ss.L; end; X0 = [K0, L0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, GDense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: (r = %.4f, K = %.2f),  markets(S = %.2d,  Y-C-I = %.2d, Kgap = %.2d, Lgap = %.2d) \n\n', ...
        ss.r, ss.K, ss.S, ss.excess_supply, ss.excess_capital, ss.excess_labor);
    
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
    
    figure('visible', 'off'); hold on;
    l1=scatter(G_adapt{n}.k, V_adapt{n}(:, 1)); 
    l2=scatter(G_adapt{n}.k, V_adapt{n}(:, 2)); 
    hold off; xlabel('Capital');
    legend([l1,l2], {'$V^U(k)$', '$V^E(k)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off




