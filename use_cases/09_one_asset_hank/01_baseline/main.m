%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium and transition dynamics of
% a one-asset HANK model similar to McKay-Nakamura-Steinsson (AER, 2016) 
% using adaptive sparse grids.
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

param = define_parameters('shock_type', 'demand', 'implicit_g', 1, 'T', 30, 'N', 40, 'bfun_type', "cheb", 'cheb_H', 15);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(param.l_dense, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'a'});
G_dense.dx = G_dense.da;

% Sparse grid:
G = setup_grid(param.l, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'a'});

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
    r0 = 0.01; Y0 = 1; if exist('ss', 'var'), r0 = ss.r; Y0 = ss.Y; end; X0 = [r0, Y0]; J0 = [];
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, GDense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: (r=%.4f, Y=%.2f),  markets(S=%.2d  Y-C=%.2d  Bgap=%.2d  Lgap=%.2d) \n\n', ...
        ss.r, ss.Y, ss.S, ss.excess_supply, ss.excess_bonds, ss.excess_labor);
    
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

% MIT shock: 
switch param.shock_type
    case 'productivity'
        shock_lev = 0.05; shock_per = log(2);
    case 'monetary'
        shock_lev = 0.01; shock_per = log(2); 
    case 'demand'
        shock_lev = 0.005; shock_per = log(2); 
end
shock_t = shock_lev * ones(param.N, 1);
for n = 1:param.N-1
    shock_t(n+1) = exp(-shock_per * param.t(n+1))*shock_lev;
end

fprintf('Impulse response paths:  %.i quarters,  %.i time steps,  using %.i %s BFs\n\n', ...
         param.T, param.N, param.H(1), param.bfun_type);

% Initialize paths and grid: (guessing paths for Y and L)
X0 = [ss.N, ss.r] .* ones(param.N, 1);
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
    l1=scatter(G_adapt{n}.a, V_adapt{n}(:, 1)); 
    l2=scatter(G_adapt{n}.a, V_adapt{n}(:, 2)); 
    hold off; xlabel('Wealth');
    legend([l1, l2], {'$V^U(a)$', '$V^E(a)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');    
    set(gcf, 'renderer', 'Painters');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

% Transition:
figure('visible', 'off');
subplot(3, 3, 1); 
plot(sim.t, (sim.Y - ss.Y)/ss.Y); ylabel('% dev'); title('$Y_t$', 'Interpreter', 'Latex');
subplot(3, 3, 2); 
plot(sim.t, (sim.N - ss.N)/ss.N); title('$N_t$', 'Interpreter', 'Latex');
subplot(3, 3, 3); 
plot(sim.t, (sim.w - ss.w)/ss.w); title('$w_t$', 'Interpreter', 'Latex');
subplot(3, 3, 4); 
plot(sim.t, sim.i); ylabel('lvl'); title('$i_t$', 'Interpreter', 'Latex');
subplot(3, 3, 5); 
plot(sim.t, sim.pi); title('$\pi_t$', 'Interpreter', 'Latex');
subplot(3, 3, 6); 
plot(sim.t, sim.r); title('$r_t$', 'Interpreter', 'Latex');
subplot(3, 3, 7); 
plot(sim.t, (sim.Pi - ss.Pi)/ss.Pi); ylabel('% dev'); xlabel('Quarters'); title('$\Pi_t$', 'Interpreter', 'Latex');
subplot(3, 3, 8); 
plot(sim.t, (sim.tau - ss.tau)/ss.tau); xlabel('Quarters'); title('$\tau_t$', 'Interpreter', 'Latex');
subplot(3, 3, 9); 
plot(sim.t, shock_t); xlabel('Quarters'); title([param.shock_type, ' shock (lvl)'], 'Interpreter', 'Latex');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/transition_dynamics.eps');


diary off




