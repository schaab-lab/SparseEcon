function g = KF(s_dense, G, param)


%% RECONSTRUCT OPERATORS
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

At = FD_operator(G, ones(G.J, 1), zeros(G.J, 1), 2, '1');

Aa{1} = FD_operator(G, s_dense{1}, zeros(G.J, 1), 1, '1');
Aa{2} = FD_operator(G, s_dense{2}, zeros(G.J, 1), 1, '2');

AT = (blkdiag(Aa{1} + At, Aa{2} + At) + Az)';


%% POPULATION DYNAMICS
% g = zeros(G.J, param.discrete_types);
% g(G.a == param.amin, :) = 1/numel(param.zz) / G.dx / sum(G.a==param.amin);
% g = ones(G.J, param.discrete_types);
% g = g / (numel(g) * G.dx);
g = zeros(G.J, param.discrete_types);
g(G.grid(:, 1)<0.2, :) = 1/ ( numel(g(G.grid(:, 1)<0.2, :)) * G.dx);


%{
% Initialize births at 0 wealth:
[~, birth_idx] = min( (G.a).^2 + (G.t-param.tmin).^2);
if G.a(birth_idx, :) < 0, birth_idx2 = birth_idx+1; elseif G.a(birth_idx, :) > 0, birth_idx2 = birth_idx-1; end
birth_ID = zeros(G.J, 2); 
if G.a(birth_idx) == 0 
    birth_ID(birth_idx, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)];
else
    birth_ID(birth_idx, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)] ...
                            * ( abs(G.a(birth_idx2)) / ...
                              ( abs(G.a(birth_idx)) + abs(G.a(birth_idx2))));
    birth_ID(birth_idx2, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)] ...
                            * ( abs(G.a(birth_idx)) / ...
                              ( abs(G.a(birth_idx)) + abs(G.a(birth_idx2))));
end

% Mass of newborn cohort:
% conceptually: kappa = sum(sum( g(G.t == param.tmax) * G.dx ));
% kappa = 0; %sum(sum( g(G.t == param.tmax) * G.dx ));
kappa = - sum(AT * [g(:, 1); g(:, 2)]); 
%}
% param.Delta_KF = 0.5;
% param.maxit_KF = 10000;

%% SOLVE KF PDE
for n = 1:param.maxit_KF
    
    % EXPLICIT:
    % 1/D *(gnew - g) = AT * g + kappa * birth_ID ;
    % gnew = g + D * AT * g + D * kappa * birth_ID;    
    % g_new = [g(:, 1); g(:, 2)] + param.Delta_KF * AT * [g(:, 1); g(:, 2)]; % + param.Delta_KF * kappa * [birth_ID(:, 1); birth_ID(:, 2)];
    
    % IMPLICIT:
    % 1/D *(gnew - g) = AT * gnew + kappa * birth_ID
    % (1/D - AT) * gnew = 1/D*g + kappa*birth_ID
    B = 1/param.Delta_KF .* speye(param.discrete_types*G.J) - AT;
    b = [g(:, 1); g(:, 2)] / param.Delta_KF; % + kappa * [birth_ID(:, 1); birth_ID(:, 2)];

    g_new = B\b;

    diff = max(abs( [g(:, 1); g(:, 2)] - g_new ));
    if diff < param.crit_KF, break; end
    g = [g_new(1:G.J), g_new(1+G.J:end)];    
end
if n == param.maxit_KF, fprintf('KF did not converge. Remaining Gap: %.2d\n', diff); end


% Some tests:
mass = sum(g * G.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end


end