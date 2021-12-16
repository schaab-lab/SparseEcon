%------------------------------------------------------------------------%
% 
% This baseline extends the code written by Greg Kaplan and Ben Moll (see
% https://benjaminmoll.com) based on Skiba (1978) to adaptive sparse grids.
% 
% Code written by Andreas Schaab and Allen T. Zhang.
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

param = define_parameters();


%% STEADY STATE CAPITAL
kH = (param.alpha * param.AH / (param.rho+param.delta))^(1/(1-param.alpha)) + param.kappa;
Kss = param.kappa / (1 - (param.AL/param.AH)^(1/param.alpha));

param.kmin = 0.01 * kH;
param.kmax = 1.30 * kH;

param.min(1) = param.kmin;
param.max(1) = param.kmax;


%% INITIALIZE GRIDS
G = setup_grid(param.l, param.surplus, param.min, param.max, ...
    'NamedDims', {1}, 'Names', {'k'});


%% PRODUCTION FUNCTION
G.y = skiba_production_function(G, param);


%% DYNAMIC PROGRAMMING PROBLEM
blacklist = []; J0 = [];
V_adapt = cell(param.max_adapt_iter, 1); 
s_adapt = cell(param.max_adapt_iter, 1); 
c_adapt = cell(param.max_adapt_iter, 1); 
G_adapt = cell(param.max_adapt_iter, 1); 

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE VFI    
    G.income = G.y - param.delta * G.k;

    % State-constrained boundary conditions:
    left_bound  = param.u1(G.income);
    right_bound = param.u1(G.income);

    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) sparse_project(left_bound,  points, G);
    BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
    G = gen_FD(G, BC);

    % Initialize guess V0:
    if ~isfield(G, 'V0'), G.V0 = param.u(G.income) / param.rho; end

    % Solve VFI:
    [V, hjb] = VFI(G, param);
    G.V0 = V;
    
    V_adapt{adapt_iter} = V; G_adapt{adapt_iter} = G;
    c_adapt{adapt_iter} = hjb.c; s_adapt{adapt_iter} = hjb.s;
    
    
    %% ADAPT GRID
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change == 0, break; end
    
    % Update grid objects:
    G.V0 = BH_adapt * G.V0;
    G.y  = skiba_production_function(G, param);

    
end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

% Production function:
figure('visible', 'off');
plot(G.k, G.y); xlabel('$k$', 'Interpreter', 'Latex'); xlabel('$f(k)$', 'Interpreter', 'Latex');
exportgraphics(gcf, './output/skiba_production_function.eps');

% Value function:
for n = 1:adapt_iter
    
    figure('visible', 'off'); 
    
    subplot(1, 2, 1); hold on;
    yyaxis left;  l1 = scatter(G_adapt{n}.k, V_adapt{n}); 
    yyaxis right; l2 = scatter(G_adapt{n}.k, s_adapt{n});
    hold off; xlabel('Capital: $k$', 'Interpreter', 'Latex');
    legend([l1, l2], {'$V(k)$', '$s(k)$'}, 'box', 'off', 'Interpreter', 'Latex', 'Location', 'NorthWest');    
    
    subplot(1, 2, 2); hold on;
    l1 = scatter(G_adapt{n}.k, c_adapt{n});
    l2 = scatter(G_adapt{n}.k, G_adapt{n}.y - param.delta * G_adapt{n}.k);
    hold off; xlabel('Capital: $k$', 'Interpreter', 'Latex');
    legend([l1, l2], {'$c(k)$', '$f(k) - \delta k$'}, 'box', 'off', 'Interpreter', 'Latex', 'Location', 'NorthWest');    

    set(gcf, 'renderer', 'Painters');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off




