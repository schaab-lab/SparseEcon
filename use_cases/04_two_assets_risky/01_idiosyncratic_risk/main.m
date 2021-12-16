%------------------------------------------------------------------------%
% 
% This code solves a two-asset version of the Aiyagari (1994) model using
% adaptive sparse grids. One asset is liquid and riskfree (bonds), the 
% other is liquid but risky (stocks). Asset return risk is idiosyncratic. 
% We build on code written by Ben Moll (see https://benjaminmoll.com).
% 
% Code written by Andreas Schaab and Allen T. Zhang.
% Current version: December 2021. First version: February 2020.
% 
% If you find this code helpful in your own work, please cite:
%   Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%   with Adaptive Sparse Grids. Working Paper.
% Thanks!
%
% Replication output:
% (r=0.0151  K=6.75)  markets(S=4.39e-08  B=-3.78e-14  Y-C-I=-2.79e-13  K-KS=-1.97e-06)
% [r=0.0151  K=6.75]  Markets: [S=4.39e-08  B=-2.54e-11  Y-C-I=-7.19e-12  K-KS=-1.97e-06]
% 
% Switching to asset-pricing BC at nmax (away from Huggett), so saving possible at nmax!
% [r=0.0151  K=6.75]  Markets: [S=4.39e-08  B=-5.47e-14  Y-C-I=1.71e-13  K-KS=-1.97e-06] 
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
G_dense = setup_grid(param.l_dense, 0, param.min, param.max, ...
    'NamedDims', {1}, 'Names', {'n'}, 'DxxDims', 1);
G_dense.dx = G_dense.dn;

% Sparse grid:
G = setup_grid(param.l, 0, param.min, param.max, ...
    'NamedDims', {1}, 'Names', {'n'}, 'DxxDims', 1);

% Initialize boundary conditions for each discrete type:
for j = 1:param.discrete_types
    BC{1}.left.type = '0'; BC{1}.right.type = '0';
    G = gen_FD_chain_rule(G, BC, num2str(j));
    G_dense = gen_FD_chain_rule(G_dense, BC, num2str(j));
end

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.04; K0 = 4; if exist('ss','var'), r0 = ss.r; K0 = ss.K; end
    X0 = [r0, K0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);

    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0,[numel(X0),1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);

    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);

    fprintf('Stationary Equilibrium: [r=%.4f  K=%.2f]  Markets: [S=%.2d  B=%.2d  Y-C-I=%.2d  K-KS=%.2d] \n\n', ...
        ss.r, ss.K, ss.excess_savings, ss.excess_bonds, ss.excess_supply, ss.excess_capital);
    
    V_adapt{adapt_iter} = ss.V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change==0, break; end
    
    % Update grid objects:
    % G.V0 = BH_adapt * G.V0;
    G = rmfield(G, 'V0');
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
    
    % Update BCs:
    for j = 1:param.discrete_types
        BC{1}.left.type = '0'; BC{1}.right.type = '0';
        G = gen_FD_chain_rule(G, BC, num2str(j));
    end
    
end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n',run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off'); hold on;
    l1=scatter(G_adapt{n}.n, V_adapt{n}(:,1)); 
    l2=scatter(G_adapt{n}.n, V_adapt{n}(:,2)); 
    hold off; xlabel('Net worth');
    legend([l1,l2],{'$V^U(n)$','$V^E(n)$'},'Interpreter','Latex','box','off','Location','SouthEast');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off




