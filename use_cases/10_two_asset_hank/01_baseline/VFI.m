function [V, hjb] = VFI(G, param)

V = G.V0;

% Exogenous operators:
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for iter = 1:param.maxit

% Policy functions:
hjb = HJB(V, G, param);
if any(any(isnan(hjb.c))), V = NaN(1); return; end

% Finite-difference operators:
Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
for j = 1:param.discrete_types
    Asc{j} = FD_operator(G, hjb.sc(:, j),   zeros(G.J, 1), 1, num2str(j));
    Asi{j} = FD_operator(G, hjb.si(:, j),   zeros(G.J, 1), 1, num2str(j));
    Ami{j} = FD_operator(G, hjb.iota(:, j), zeros(G.J, 1), 2, num2str(j));
    Amk{j} = FD_operator(G, G.income_k,     zeros(G.J, 1), 2, num2str(j));
end

A = blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az;

B = (1/param.Delta + param.rho + param.deathrate)*speye(2*G.J) - A;
b = hjb.u(:) + V(:) / param.Delta;

% Solve linear system:
V_new = B \ b;
% [V_new, flag] = gmres(B, b, [], param.crit/10, 200, [], [], V(:));

% Update:
V_change = V_new - V(:);
V = [V_new(1:G.J), V_new(1+G.J:end)];

dist = max(max(abs(V_change)));
if dist < param.crit, break; end

% if mod(iter, 1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end
if ~isreal(V), fprintf('Complex values in VFI: terminating process.'); V = NaN(1); return; end

end

hjb.A = A;
if iter == param.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
