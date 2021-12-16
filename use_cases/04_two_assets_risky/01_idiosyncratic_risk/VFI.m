function [V, hjb] = VFI(G, agg, param)

V = G.V0;

% Exogenous Operators
Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];


for iter = 1:param.maxit

% COMPUTE POLICY FUNCTIONS
hjb = HJB(V, G, param);

% ASSEMBLE FD OPERATOR MATRIX
sig = param.sigR * hjb.theta .* G.n;

[An1, BC_const1] = FD_operator(G, hjb.s(:, 1), sig(:, 1), 1, '1');
[An2, BC_const2] = FD_operator(G, hjb.s(:, 2), sig(:, 2), 1, '2');

assert(all(all( hjb.s(G.n == param.nmin,:) >= 0 ))); %=> VnB BC doesn't mess with A
% We no longer use Huggett BC at nmax, so household can save at nmax!
% assert(all(all( hjb.s(G.n == param.nmax,:) <= 0 ))); %=> VnF BC doesn't mess with A
assert(all(all( sig(G.n == param.nmin,:)   == 0 ))); %=> Vnn BC doesn't mess with A at nmin

A = blkdiag(An1, An2) + Az; BC_const = [BC_const1; BC_const2];

B = (1/param.Delta + param.rho)*speye(2*G.J) - A;
b = [hjb.u(:, 1); hjb.u(:, 2)] + [V(:, 1); V(:, 2)] / param.Delta + BC_const;

% SOLVE LINEAR SYSTEM
V_new = B\b;
% [V_new,flag] = gmres(B, b, [], param.crit/10, 500, [], [], [V(:,1); V(:,2)]);

% UPDATE
V_change = V_new - [V(:,1); V(:,2)];
V = [V_new(1:G.J), V_new(1+G.J:end)];

dist = max(max(abs(V_change)));
if dist < param.crit, break; end

% if mod(iter,1)==0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end
if ~isreal(V), fprintf('Complex values in VFI: terminating process.'); V = NaN(1); return; end

end

hjb.A = A;
if iter == param.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
