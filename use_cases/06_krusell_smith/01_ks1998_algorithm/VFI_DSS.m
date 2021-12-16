function [V, hjb] = VFI_DSS(G, param)

V = G.V0;

% EXOGENOUS OPERATORS
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for iter = 1:param.maxit

% COMPUTE POLICY FUNCTIONS
hjb = HJB(V, G, param);
if any(any(isnan(hjb.c))), V = NaN(1); return; end

% ASSEMBLE FD OPERATOR MATRIX
Aa{1} = FD_operator(G, hjb.s(:, 1), zeros(G.J, 1), 1, '1');
Aa{2} = FD_operator(G, hjb.s(:, 2), zeros(G.J, 1), 1, '2');

A = blkdiag(Aa{1}, Aa{2}) + Az;

B = (1/param.Delta + param.rho)*speye(2*G.J) - A;
b = [hjb.u(:, 1); hjb.u(:, 2)] + [V(:, 1); V(:, 2)] / param.Delta;

% SOLVE LINEAR SYSTEM
V_new = B\b;
% [V_new, flag] = gmres(B, b, [], param.crit/10, 2000, [], [], [V(:, 1); V(:, 2)]);

% UPDATE
V_change = V_new - [V(:, 1); V(:, 2)];
V = [V_new(1:G.J), V_new(1+G.J:end)];

dist = max(max(abs(V_change)));
if dist < param.crit, break; end

% if mod(iter,1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end
if ~isreal(V), fprintf('Complex values in VFI: terminating process.'); V = NaN(1); return; end

end

hjb.A = A;
if iter == param.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
