function [V, hjb] = VFI(G, agg, param)

V = G.V0;

% Exogenous Operators
% Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
%        speye(G.J)*param.la2, -speye(G.J)*param.la2];
[Az, const_z] = FD_operator(G, (param.theta_z*(param.zmean - log(G.z)) + param.sig_z^2/2).*G.z, param.sig_z*G.z, 2);


for iter = 1:param.maxit

% COMPUTE POLICY FUNCTIONS
hjb = HJB(V, G, param);
if any(any(isnan(hjb.c))), V = NaN(1); return; end

% ASSEMBLE FD OPERATOR MATRIX
% AA = [];
% for j = 1:param.discrete_types
%     A = FD_operator(G, hjb.mu{j}(:,1), hjb.sig{j}(:,1), 1);
%     for i = 2:param.d_idio
%         A = A + FD_operator(G, hjb.mu{j}(:,i), hjb.sig{j}(:,i), i);
%     end
%     AA = blkdiag(AA, A);
% end
% A = AA + Az;
[Aa, const_a] = FD_operator(G, hjb.s, zeros(G.J,1), 1);

A = Aa + Az;
const = const_a + const_z;

B = (1/param.Delta + param.rho)*speye(G.J) - A;
b = hjb.u + V / param.Delta + const;

% SOLVE LINEAR SYSTEM
V_new = B\b;
% [V_new, flag] = gmres(B, b, [], param.crit/10, 2000, [], [], [V(:,1); V(:,2)]);

% UPDATE
V_change = V_new - V;
V = V_new;

dist = max(max(abs(V_change)));
if dist < param.crit
    break
end

% if mod(iter,1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end

if ~isreal(V), fprintf('Complex values in VFI: terminating process.'); V = NaN(1); return; end

end

hjb.A = A;
if iter == param.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
