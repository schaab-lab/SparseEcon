%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium and transition dynamics in 
% response to a productivity shock for a model similar to Aiyagari (1994)
% using adaptive sparse grids.
% 
% Code written by Andreas Schaab and Allen Zhang.
% Current version: September 2021. First version: December 2019.
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
warning off

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
    'NamedDims', {1}, 'Names', {'k'});
G_dense.dx = G_dense.dk;

% Sparse grid:
G = setup_grid(param.l, 0, param.min, param.max, ...
    'NamedDims', {1}, 'Names', {'k'});

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
fprintf('\n\n:::::::::::   STATIONARY EQUILIBRIUM   ::::::::::: ');

blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  Grid Adaptation Iteration %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    K0 = 22; L0 = 1; if exist('ss', 'var'), K0 = ss.K; L0 = ss.L; end; X0 = [K0, L0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: (r=%.4f, K=%.2f),  markets(S=%.2d  Y-C-I=%.2d  Kgap=%.2d  Lgap=%.2d) \n\n', ...
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


%% COMPUTE TRANSITION DYNAMICS
fprintf('\n\n:::::::::::   TRANSITION DYNAMICS   ::::::::::: \n\n');

% Productivity shock: 
shock_lev = 0.05;
shock_per = log(2);
shock_t = shock_lev * ones(param.N, 1);
for n = 1:param.N-1
    shock_t(n+1) = exp(-shock_per * param.t(n+1))*shock_lev;
end

fprintf('Impulse response paths:  %.i quarters,  %.i time steps,  using %.i %s BFs\n\n', ...
         param.T, param.N, param.H(1), param.bfun_type);

% Initialize paths and grid: (guessing paths for Y and L)
X0 = [ss.Y, ss.L] .* ones(param.N, 1);
[PHI0, param.nodes] = basis_fun_irf(X0, [], param.H(1), param.H(2), ...
    param.bfun_type, param.t, "get_coefficient");

[diff0, G, G_dense, ~] = transition(PHI0, G, G_dense, shock_t, ss, param);

% Solve for prices:
f = @(x, y) transition(x, y{1}, y{2}, shock_t, ss, param); y0{1} = G; y0{2} = G_dense;
PHI = fsolve_newton(f, reshape(PHI0, [numel(PHI0), 1]), diff0, y0, 0, 5, 2);

% Update everything given new prices:
[diff, G, G_dense, sim] = transition(PHI, G, G_dense, shock_t, ss, param);
sim.PHI = PHI; sim.param = param;


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

% DSS: 
fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off'); hold on;
    l1=scatter(G_adapt{n}.k, V_adapt{n}(:, 1)); 
    l2=scatter(G_adapt{n}.k, V_adapt{n}(:, 2)); 
    hold off; xlabel('Capital');
    legend([l1, l2], {'$V^U(k)$', '$V^E(k)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

% Transition:
figure('visible', 'off');
subplot(2, 2, 1); 
plot(sim.t, (sim.Y - ss.Y)/ss.Y); ylabel('% dev'); title('$Y_t$', 'Interpreter', 'Latex');
subplot(2, 2, 2); 
plot(sim.t, (sim.K - ss.K)/ss.K); title('$K_t$', 'Interpreter', 'Latex');
subplot(2, 2, 3); 
plot(sim.t, (sim.L - ss.L)/ss.L); ylabel('% dev'); xlabel('Quarters'); title('$L_t$', 'Interpreter', 'Latex');
subplot(2, 2, 4); 
plot(sim.t, exp(sim.Z)); xlabel('Quarters'); title('$Z_t$', 'Interpreter', 'Latex');
exportgraphics(gcf, './output/transition_dynamics.eps');



diary off




