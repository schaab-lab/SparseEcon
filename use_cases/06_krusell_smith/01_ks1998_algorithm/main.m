%------------------------------------------------------------------------%
% 
% This code solves the seminal Krusell-Smith (1998) model using adaptive 
% sparse grids. The code uses the original Krusell-Smith algorithm and
% produces a Den Haan (2010) metric of 0.14%.
% 
% Code written by Andreas Schaab.
% Current version: December 2021. First version: September 2019.
% 
% If you find this code helpful in your own work, please cite:
%  1. Schaab, A. Micro and Macro Uncertainty. Working Paper.
%  2. Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time 
%     with Adaptive Sparse Grids. Working Paper.
% Thanks!
% 
% Replication check:
%  Stationary equilibrium
%     (r = 0.0495,  K = 6.13)   markets(S = 1.58e-09,  Y-C-I = 1.58e-09, Kgap = 4.05e-13) 
%  Final print after 65 KS iterations
%     Den Haan (2010) metric for KS algorithm (max): 0.142828 
%     Den Haan (2010) metric for KS algorithm (avg): 0.021446 
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

param = define_parameters('max_adapt_iter', 1, 'T', 2000, 'N', 5000, 'lambda_LOM', 0.85);


%% INITIALIZE GRIDS

% Dense grid:
G_dense = setup_grid(param.l_dense, 0, param.min(1), param.max(1), ...
    'NamedDims', {1}, 'Names', {'a'});
G_dense.dx = G_dense.da;

% Sparse grid:
% G = setup_grid(param.l, param.surplus, param.min(1), param.max(1), ...
%     'NamedDims', {1}, 'Names', {'a'});
G = G_dense;

% Initialize boundary conditions:
% BC{1}.left.type = '1'; BC{1}.right.type = '1';
% G_dense = gen_FD(G_dense, BC);
% G = gen_FD(G, BC);

% Projection matrix:
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);


%% COMPUTE (D)ETERMINISTIC (S)TEADY (S)TATE ON ADAPTED SPARSE GRID
blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1); 
c_adapt = cell(param.max_adapt_iter, 1); 
s_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE STATIONARY EQUILIBRIUM
    r0 = 0.04; if exist('ss', 'var'), r0 = ss.r; end; X0 = r0; J0 = [];
    
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
    
    V_adapt{adapt_iter} = ss.V; c_adapt{adapt_iter} = ss.c; s_adapt{adapt_iter} = ss.s; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    if adapt_iter == param.max_adapt_iter, break; end
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, ss.V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
    
end


%% UPDATE STATE SPACE
param.Kmax = 1.10 * ss.K;
param.Kmin = 0.90 * ss.K;
param.K0   = ss.K;

param.max(3) = param.Kmax;
param.min(3) = param.Kmin;


%% HOUSEHOLD GRID
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1, 2, 3, 2:3}, 'Names', {'a', 'Z', 'K', 'X'}, 'DxxDims', 2);
G = update_grid(G, G.grid(G.lvl(:, 1) <= param.l_dense, :), G.lvl(G.lvl(:, 1) <= param.l_dense, :));

% Projection matrices:
% G.BH_dense  = H_basis(G_dense.grid, G_dense.lvl, G.grid(:, 1), G.lvl(:, 1)) * G.H_comp;
% G.BH_interp = agg_interp(G.grid, G.lvl, G_dense.grid, G_dense.lvl);
% G.BH_comp   = G.BH_interp * G.H_comp;


%% AGGREGATE GRID
[agg, agg.full2agg, agg.agg2full] = keep_dims(G, param.d_idio+1:param.d);

% Boundary conditions:
BC = cell(agg.d, 1);
for k = 1:agg.d
    BC{k}.left.type = '0'; BC{k}.right.type = '0';
end
agg = gen_FD(agg, BC);

% Macroeconomic aggregates:
agg.muZ  = param.thetaZ * (param.Zmean - agg.Z);
agg.sigZ = param.sigmaZ * ones(agg.J, 1);
agg.muK  = zeros(agg.J, 1);
agg.sigK = zeros(agg.J, 1);
agg.muX  = [agg.muZ, agg.muK];
agg.sigX = [agg.sigZ, agg.sigK];

agg.r = param.alpha     .* exp(agg.Z) .* agg.K .^(param.alpha-1) .* param.L.^(1-param.alpha) - param.delta;
agg.w = (1-param.alpha) .* exp(agg.Z) .* agg.K .^(param.alpha)   .* param.L.^( -param.alpha);
agg.Y = exp(agg.Z) .* agg.K.^param.alpha .* param.L.^(1-param.alpha);

if min(agg.r)<0, fprintf('New grid implies negative interest rate.'); return; end


%% KS SIM-EST LOOP

% Boundary conditions: (they stay fixed)
G.income = agg.r(agg.agg2full) .* G.a + agg.w(agg.agg2full) .* param.zz;

clear BC;
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);
for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) sparse_project(left_bound(:, j),  points, G);
    BC{1}.right.f = @(points) sparse_project(right_bound(:, j), points, G);
    for k = 2:param.d
        BC{k}.left.type = '0'; BC{k}.right.type = '0';
    end
    G = gen_FD(G, BC, num2str(j));
end

% Initialize VF
v0 = ss.V(ismember(G_dense.a, unique(G.a)), :);
[~, mapping] = ismember(G.a, unique(G.a));
G.V0 = v0(mapping, :);

for outer = 1:param.max_KS
    fprintf('\n\n ---------   KS ITERATION #%.i   ---------- \n', outer);
    
    % VFI
    [V, hjb] = VFI(G, agg, param);
    
    % Simulate
    sim = sim_fun(G, G_dense, agg, V, ss, param);
    
    % Estimate
    [basis_sim, basis_agg] = estimation_model(sim, agg, param.estimation_model_type, 1:2);     
    dK = diff(sim.K(sim.n_data)); %./ sim.K(sim.n_data(1:end-1));           
    PHI_muK = regress(dK, basis_sim(sim.n_data(1:end-1), :) * sim.dt);
    new_muK = basis_agg * PHI_muK;
    
    fprintf('Remaining Update Dist = %.2d \n', max(abs(agg.muK-new_muK)));
    fprintf('Remaining DH Metric   = %.3f \n', 100 * max(abs( log(sim.K) - log(sim.K_lom) )));
    if max(abs(agg.muK-new_muK)) < param.crit_KS, break; end
    
    % Update
    agg.muK = param.lambda_LOM * agg.muK + (1-param.lambda_LOM) * new_muK;
    
    if nnz(sim.GK == 0 | sim.GK == 1) / numel(sim.GK) > 0
        fprintf('\nProportion of K on boundary: KF  simulation %g.\n', ...
            nnz(sim.GK == 0 | sim.GK == 1) / numel(sim.GK));
    elseif nnz(sim.GK_lom == 0 | sim.GK_lom == 1) / numel(sim.GK_lom) > 0
        fprintf('\nProportion of K on boundary: LOM simulation %g.\n', ...
            nnz(sim.GK_lom == 0 | sim.GK_lom == 1) / numel(sim.GK_lom));
    end
    
end

fprintf('Hash #1: %.7f,  Hash #2: %.7f \n\n', mean(sim.K), mean(sim.K_lom));


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

% DSS:
fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off'); hold on;
    l1 = scatter(G_adapt{n}.a, V_adapt{n}(:, 1)); 
    l2 = scatter(G_adapt{n}.a, V_adapt{n}(:, 2)); 
    hold off; xlabel('Capital');
    legend([l1, l2], {'$V^U(a)$', '$V^E(a)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');
    set(gcf, 'renderer', 'Painter');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end


% Simulation:
lStyle  = {'-', '-.'};
figure('visible', 'off'); hold on;
l1 = plot(param.t, sim.K,     'LineStyle', lStyle{1}, 'Color', colorPalette(1, :), 'LineWidth', 2.5);
l2 = plot(param.t, sim.K_lom, 'LineStyle', lStyle{2}, 'Color', colorPalette(2, :), 'LineWidth', 2);
hold off; xlabel('Simulation horizon (years)'); ylabel('Capital');
set(gcf, 'units', 'inches', 'position', [5.5, 4.5, 10, 5.5])
legend([l1, l2], {'Capital: simulation', 'Capital: forecast'}, ...
    'Location', 'SouthOutside', 'NumColumns', 2, 'FontSize', 12, 'Interpreter', 'Latex', 'box', 'off');
set(gcf, 'renderer', 'Painter');
exportgraphics(gcf, './output/simulation_dh.eps');


% Den Hann (2010) metric:
denhaan2010 = 100 * max(abs( log(sim.K) - log(sim.K_lom) ));
denhaan2010avg = 100 * mean(abs( log(sim.K) - log(sim.K_lom) ));
fprintf('The Den Haan (2010) metric for KS algorithm (max): %.6f \n', denhaan2010);
fprintf('The Den Haan (2010) metric for KS algorithm (avg): %.6f \n', denhaan2010avg);


diary off

