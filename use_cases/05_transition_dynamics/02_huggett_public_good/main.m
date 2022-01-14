%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium and transition dynamics in 
% a variant of the Huggett model with diffusive earnings risk and a public
% good.
% 
% Code written by Andreas Schaab.
% Current version: January 2022. First version: January 2022.
% 
% If you find this code helpful in your own work, please cite:
% - Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
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

param = define_parameters('max_adapt_iter', 1, 'T', 20, 'N', 60);


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
fprintf('\n\n:::::::::::   STATIONARY EQUILIBRIUM   ::::::::::: ');

blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  Grid Adaptation Iteration %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.0018; if exist('ss', 'var'), r0 = ss.r; end; X0 = r0;
    
    % Get better guess for value function:
    [diff0, G, G_dense, ~] = stationary(X0, G, G_dense, param);
    
    % Solve for steady state prices:
    % f = @(x, y) stationary(x, y, G_dense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', false, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x, G, G_dense, param), X0, options);
    
    % Solve with correct prices:
    [~, G, G_dense, ss] = stationary(X, G, G_dense, param);
    
    fprintf('Stationary Equilibrium: r = %.4f,  markets(B = %.2d,  S = %.2d,  Y-C-G = %.2d) \n\n', ...
        ss.r, ss.B, ss.S, ss.excess_supply);
    
    V_adapt{adapt_iter} = ss.V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    if adapt_iter == param.max_adapt_iter, break; end
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
shock = param.shock_level * ones(param.N, 1);
for n = 1:param.N-1
    shock(n+1) = exp(-param.theta_shock * param.t(n+1))*param.shock_level;
end

% Policy shock: 
policy = param.policy + param.policy_shock * ones(param.N, 1);
for n = 1:param.N-1
    policy(n+1) = param.policy + exp(-param.theta_policy * param.t(n+1))*param.policy_shock;
end

fprintf('Impulse response paths:  %.i quarters,  %.i time steps,  using %.i %s BFs\n\n', ...
         param.T, param.N, param.H(1), param.bfun_type);

% Initialize paths and grid: (guessing path for r)
X0 = ss.r .* ones(param.N, 1);
[PHI0, param.nodes] = basis_fun_irf(X0, [], param.H(1), param.H(2), ...
    param.bfun_type, param.t, "get_coefficient");

[diff0, G, G_dense, ~] = transition(PHI0, G, G_dense, shock, policy, ss, param);

% Solve for prices:
f = @(x, y) transition(x, y{1}, y{2}, shock, policy, ss, param); y0{1} = G; y0{2} = G_dense;
PHI = fsolve_newton(f, reshape(PHI0, [numel(PHI0), 1]), diff0, y0, 0, 5, 2);

% Update everything given new prices:
[diff, G, G_dense, sim] = transition(PHI, G, G_dense, shock, policy, ss, param);
sim.PHI = PHI; sim.param = param;


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

% DSS: 
for n = 1:adapt_iter
    
    figure('visible', 'off');
    l1 = scatter3(G_adapt{n}.a, G_adapt{n}.z, V_adapt{n}(:, 1));
    xlabel('Wealth: $a$', 'Interpreter', 'Latex');
        xlh = get(gca, 'xlabel'); gxl = get(xlh); xlp = get(xlh, 'Position');
        set(xlh, 'Rotation', 16, 'Position', xlp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    ylabel('Earnings: $z$', 'Interpreter', 'Latex');
        ylh = get(gca, 'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh,'Rotation', -27, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    zlabel('$V(a,z)$', 'Interpreter', 'Latex');
        zlh = get(gca, 'zlabel'); gzl = get(zlh); zlp = get(zlh, 'Position');
        set(zlh, 'Rotation', 0, 'Position', zlp+[0, 0.08, 0], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    set(gcf, 'renderer', 'Painters');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

% Transition dynamics:
figure('visible', 'off');
subplot(2, 2, 1);
plot(sim.t, 100 * (sim.Y - ss.Y)/ss.Y); ylabel('% dev'); title('$Y_t$', 'Interpreter', 'Latex');
subplot(2, 2, 2); 
plot(sim.t, 100 * (sim.G - param.policy)/param.policy); title('$G_t$', 'Interpreter', 'Latex');
subplot(2, 2, 3);
plot(sim.t, sim.r); ylabel('lvl'); xlabel('Quarters'); title('$r_t$', 'Interpreter', 'Latex');
subplot(2, 2, 4);
plot(sim.t, exp(sim.Z)); xlabel('Quarters'); title('$Z_t$', 'Interpreter', 'Latex');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/transition_dynamics.eps');

% Transition market clearing:
figure('visible', 'off');
subplot(1, 3, 1);
plot(sim.t, sim.excess_bonds); ylabel('level'); title('$B_t$', 'Interpreter', 'Latex');
subplot(1, 3, 2); 
plot(sim.t, sim.excess_goods); xlabel('Quarters'); title('$Y_t - C_t - G_t$', 'Interpreter', 'Latex');
subplot(1, 3, 3); 
plot(sim.t, sim.excess_saving); title('$S_t$', 'Interpreter', 'Latex');
set(gcf, 'renderer', 'Painters');
exportgraphics(gcf, './output/transition_market_clearing.eps');


diary off




