function g = KF(sc, si, mi, mk, G, param)

% RECONSTRUCT A MATRIX:
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
for j = 1:param.discrete_types
    Asc{j} = FD_operator(G, sc{j}, zeros(G.J, 1), 1, num2str(j));
    Asi{j} = FD_operator(G, si{j}, zeros(G.J, 1), 1, num2str(j));
    Ami{j} = FD_operator(G, mi{j}, zeros(G.J, 1), 2, num2str(j));
    Amk{j} = FD_operator(G, mk{j}, zeros(G.J, 1), 2, num2str(j));
end
AT = (blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az)';

% SOLVE KF PDE:
g = zeros(G.J, 2);
g(G.k == param.kmin & G.a == param.amin, :) = 1/numel(param.zz) / G.dx;

for n = 1:param.maxit_KF   
    B = 1/param.Delta_KF * speye(2*G.J) - AT;
    b = g(:) / param.Delta_KF;

    g_new = B \ b;

    diff = max(abs( [g(:, 1); g(:, 2)] - g_new ));
    if diff < param.crit_KF
        break
    end
    g = [g_new(1:G.J), g_new(1+G.J:end)];    
end
if n == param.maxit_KF, fprintf('KF1 did not converge. Remaining Gap: %.2d\n', diff); end


end