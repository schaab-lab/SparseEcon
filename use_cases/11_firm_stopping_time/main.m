%------------------------------------------------------------------------%
% 
% This code solves the optimal stopping problem of a firm that can decide
% when to shut down a production plant.
% We build on code written by Ben Moll (see https://benjaminmoll.com).
% 
% Code written by Andreas Schaab.
% Current version: May 2022. First version: July 2021.
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
G = setup_grid(param.l, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'x'}, 'DxxDims', 1);

% Boundary conditions:
BC{1}.left.type = '0'; BC{1}.right.type = '0';
G = gen_FD(G, BC);


%% COMPUTE VALUE FUNCTION ON ADAPTIVE SPARSE GRID
blacklist = [];
V_adapt = cell(param.max_adapt_iter, 1);
G_adapt = cell(param.max_adapt_iter, 1);

% Guess: 
z0 = zeros(G.J, 1);

for adapt_iter = 1:param.max_adapt_iter
    
    fprintf('\n\n -------  GRID ADAPTATION ITERATION %i  ------- \n\n', adapt_iter);
    
    %% SOLVE OPTIMAL STOPPING PROBLEM
    
    % Construct FD matrices:
    [Ax, const_x] = FD_operator(G, param.theta_x*ones(G.J, 1), param.sig_x*G.x, 1);
    assert(all(all(const_x == 0)));
       
    % LCP solver:
    B = param.rho*speye(G.J) - Ax;
    q = -param.u(G.x) + B*param.outside_option(G.x); 
    % For LCP form, we need: z' (Bz + q) = 0
    %  - Let z denote excess value: v - scrap
    %  - Let B = rho I - A (as usual)
    %  - Then we're left with: q = -u + Bz
    
    z = LCP(B, q, [], [], z0, 1);
        
    LCP_error = abs(z'*(B*z + q));
    if LCP_error > param.crit
        fprintf('LCP was not solved. Remaining error: %.2d \n', LCP_error); 
    end
    
    V = z + param.outside_option(G.x);    
    
    V_adapt{adapt_iter} = V; G_adapt{adapt_iter} = G;
    
    
    %% ADAPT GRID
    [G, BH_adapt, blacklist, stats] = adapt_grid(G, V, blacklist, ...
        'AddRule', param.add_rule, 'AddTol', param.add_tol, 'KeepTol', param.keep_tol);
    if stats.n_change==0, break; end
    
    % Update guess:
    z0 = BH_adapt * z;
    
    % Update BCs:
    G = gen_FD(G, BC); % BCs don't change in this model
    
    
end


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');
for n = 1:adapt_iter
    
    figure('visible', 'off'); hold on;
    l1 = scatter(G_adapt{n}.x, V_adapt{n});
    l2 = scatter(G_adapt{n}.x, param.outside_option(G_adapt{n}.x)); 
    hold off; xlabel('Productivity');
    legend([l1, l2], {'$V(x)$', '$S(x)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');
    exportgraphics(gcf, ['./output/grid_adaptation', num2str(n-1), '.eps']);

end

diary off




