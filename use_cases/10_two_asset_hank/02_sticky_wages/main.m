%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium and transition dynamics of
% a two-asset HANK model similar to Kaplan-Moll-Violante (AER, 2018) using 
% adaptive sparse grids.
% 
% Code written by Andreas Schaab.
% Current version: December 2021. First version: August 2019.
% 
% If you find this code helpful in your own work, please cite:
%  1. Schaab, A. Micro and Macro Uncertainty. Working Paper.
%  2. Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%     with Adaptive Sparse Grids. Working Paper.
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

param = define_parameters('deathrate', 0, 'add_tol', 5e-5, 'keep_tol', 1e-5, 'max_adapt_iter', 1, 'kappa', 60, ...
                          'shock_type', 'monetary', 'implicit_g', 1, 'T', 150, 'N', 180, 'bfun_type', "cheb", 'cheb_H', 30);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(0, param.l_dense, param.min, param.max, ...
    'NamedDims', {1, 2}, 'Names', {'a', 'k'});
G_dense.dx = G_dense.da*G_dense.dk;

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


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
fprintf('\n\n:::::::::::   STATIONARY EQUILIBRIUM   ::::::::::: ');

blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.03; K0 = 6; L0 = 0.75; if exist('ss', 'var'), r0 = ss.r; K0 = ss.K; L0 = ss.L; end
    X0 = [r0, K0, L0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-15);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf(['Stationary Equilibrium: [r = %.4f, K=%.4f, L=%.4f, Q=%.4f]   ', ...
             'Markets=[G=%.1d, B=%.1d, L=%.1d, K=%.1d])\n\n'], ...
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
    
    % Update BCs:
    for j = 1:param.discrete_types
        BC{1}.left.type = '0'; BC{1}.right.type = '0';
        BC{2}.left.type = '0'; BC{2}.right.type = '0';
        G = gen_FD(G, BC, num2str(j));
    end
end


%% COMPUTE TRANSITION DYNAMICS
fprintf('\n\n:::::::::::   TRANSITION DYNAMICS   ::::::::::: \n\n');

% Productivity shock: 
switch param.shock_type
    case 'productivity'
        shock_lev = 0.05; shock_per = log(2);
    case 'monetary'
        shock_lev = 0.01; shock_per = log(2); 
end
shock_t = shock_lev * ones(param.N, 1);
for n = 1:param.N-1
    shock_t(n+1) = exp(-shock_per * param.t(n+1))*shock_lev;
end

fprintf('Impulse response paths:  %.i quarters,  %.i time steps,  using %.i %s BFs\n\n', ...
         param.T, param.N, param.H(1), param.bfun_type);

% Initialize paths and grid: (guessing paths for Y and L)
X0 = [ss.r, ss.K, ss.L] .* ones(param.N, 1);
param.H(2) = 3;
[PHI0, param.nodes] = basis_fun_irf(X0, [], param.H(1), param.H(2), param.bfun_type, param.t, "get_coefficient");

[diff0, G, G_dense, ~] = transition(PHI0, G, G_dense, shock_t, ss, param);

% Solve for prices:
f = @(x, y) transition(x, y{1}, y{2}, shock_t, ss, param); y0{1} = G; y0{2} = G_dense;
PHI = fsolve_newton(f, reshape(PHI0, [numel(PHI0), 1]), diff0, y0, 0, 5, 2);

% Update everything given new prices:
[diff, G, G_dense, sim] = transition(PHI, G, G_dense, shock_t, ss, param);
sim.PHI = PHI; sim.param = param;


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

% DSS: 
for n = 1:adapt_iter
    
    figure('visible', 'off');
    l1=scatter3(G_adapt{n}.a, G_adapt{n}.k, V_adapt{n}(:, 1));
    xlabel('Liquid: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Illiquid: $k$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh, 'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    zlabel('$V^U(a, k)$', 'Interpreter', 'Latex');
        zlh = get(gca, 'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
        set(zlh, 'Rotation', 0, 'Position', zlp+[0, 11, -2], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');        
    set(gcf, 'renderer', 'Painters');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end


% Transition:
figure('visible', 'off');
subplot(3, 3, 1); 
plot(sim.t, (sim.Y - ss.Y)/ss.Y); ylabel('% dev'); title('$Y_t$', 'Interpreter', 'Latex');
subplot(3, 3, 2); 
plot(sim.t, (sim.K - ss.K)/ss.K); title('$K_t$', 'Interpreter', 'Latex');
subplot(3, 3, 3); 
plot(sim.t, (sim.L - ss.L)/ss.L); title('$L_t$', 'Interpreter', 'Latex');
subplot(3, 3, 4); 
plot(sim.t, (sim.C - ss.C)/ss.C); ylabel('% dev'); title('$C_t$', 'Interpreter', 'Latex');
subplot(3, 3, 5); 
plot(sim.t, (sim.I - ss.I)/ss.I); title('$I_t$', 'Interpreter', 'Latex');
subplot(3, 3, 6); 
plot(sim.t, (sim.w - ss.w)/ss.w); title('$w_t$', 'Interpreter', 'Latex');
subplot(3, 3, 7); 
plot(sim.t, (sim.r - ss.r)/ss.r); ylabel('% dev'); xlabel('Quarters'); ylabel('% dev'); title('$r_t$', 'Interpreter', 'Latex');
subplot(3, 3, 8); 
plot(sim.t, (sim.rk - ss.rk)/ss.rk); xlabel('Quarters'); title('$r_t^k$', 'Interpreter', 'Latex');
subplot(3, 3, 9); 
plot(sim.t, shock_t); xlabel('Quarters'); title([param.shock_type, ' shock (lvl)'], 'Interpreter', 'Latex');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/transition_dynamics.eps');


diary off


