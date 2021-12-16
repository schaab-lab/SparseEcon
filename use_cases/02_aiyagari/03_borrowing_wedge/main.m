%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a model similar to
% Aiyagari (1994) using adaptive sparse grids. Labor supply is inelastic 
% and earnings risk follows a diffusion process. There is an exogenous
% wedge between interest rates on savings and borrowing.
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

param = define_parameters('max_adapt_iter', 5);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);
G_dense.dx = G_dense.da * G_dense.dz;

% Sparse grid:
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1); 
c_adapt = cell(param.max_adapt_iter, 1); 
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
    
    V_adapt{adapt_iter} = ss.V; c_adapt{adapt_iter} = ss.c; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    if adapt_iter == param.max_adapt_iter, break; end
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
    subplot(1, 2, 1);
    l1 = scatter3(G_adapt{n}.a, G_adapt{n}.z, V_adapt{n});
    xlabel('Wealth: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    title('$V(a,z)$', 'Interpreter', 'Latex');
        
    subplot(1,2,2);
    l1 = scatter3(G_adapt{n}.a, G_adapt{n}.z, c_adapt{n});
    xlabel('Wealth: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    title('$c(a,z)$', 'Interpreter', 'Latex');
    
    set(gcf, 'renderer', 'Painters', 'Position', [10 10 1200 400]);
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end


diary off




