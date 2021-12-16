function g = KF(sc, si, mi, mk, G, param)

%% RECONSTRUCT OPERATORS
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

At = FD_operator(G, ones(G.J, 1), zeros(G.J, 1), 3, '1');

Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
for j = 1:param.discrete_types
    Asc{j} = FD_operator(G, sc{j}, zeros(G.J, 1), 1, num2str(j));
    Asi{j} = FD_operator(G, si{j}, zeros(G.J, 1), 1, num2str(j));
    Ami{j} = FD_operator(G, mi{j}, zeros(G.J, 1), 2, num2str(j));
    Amk{j} = FD_operator(G, mk{j}, zeros(G.J, 1), 2, num2str(j));
end
AT = (blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1} + At, Asc{2} + Asi{2} + Ami{2} + Amk{2} + At) + Az)';


%% POPULATION DYNAMICS
% g = zeros(G.J, param.discrete_types);
% g(G.a == param.amin, :) = 1/numel(param.zz) / G.dx / sum(G.a == param.amin);
% g = ones(G.J, param.discrete_types);
% g = g / (numel(g) * G.dx);
g = zeros(G.J, param.discrete_types);
g(G.grid(:, 1) < 0.2, :) = 1/ ( numel(g(G.grid(:, 1) < 0.2, :)) * G.dx);


[~, birth_idx] = min( (G.a).^2 + (G.k).^2 );
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


%% SOLVE KF PDE
for n = 1:param.maxit_KF
    
    % EXPLICIT:
    % 1/D *(gnew - g) = AT * g + kappa * birth_ID ;
    % g_new = g + D * AT * g + D * kappa * birth_ID;    
    % g_new = [g(:, 1); g(:, 2)] + param.Delta_KF * AT * [g(:, 1); g(:, 2)]; % + param.Delta_KF * kappa * [birth_ID(:, 1); birth_ID(:, 2)];
    
    % IMPLICIT:
    % 1/D *(g_new - g) = AT * g_new + kappa * birth_ID
    % (1/D - AT) * g_new = 1/D*g + kappa*birth_ID
    B = (1/param.Delta_KF + param.deathrate) .* speye(2*G.J) - AT;
    b = [g(:, 1); g(:, 2)] / param.Delta_KF + param.deathrate * [birth_ID(:, 1); birth_ID(:, 2)] / G.dx;

    g_new = B\b;

    diff = max(abs( [g(:, 1); g(:, 2)] - g_new ));
    if diff < param.crit_KF
        break
    end
    g = [g_new(1:G.J), g_new(1+G.J:end)];    
end
if n == param.maxit_KF, fprintf('KF1 did not converge. Remaining Gap: %.2d\n', diff); end


end