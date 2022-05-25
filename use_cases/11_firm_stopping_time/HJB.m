function hjb = HJB(V, G, param)

VaF = zeros(G.J, param.discreteTypes);
VaB = zeros(G.J, param.discreteTypes);
for j = 1:param.discreteTypes
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

cF = param.u1inv(VaF);
cB = param.u1inv(VaB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

IF = (sF>0);             % BC takes care of this: (G.grid(:,1)<1)
IB = (sB<0) .* (IF==0);  % BC takes care of this: (G.grid(:,1)>0)
I0 = (1-IF-IB);

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

% COLLECT OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

% for j = 1:param.discreteTypes, hjb.mu{j}  = s(:,j); end
% for j = 1:param.discreteTypes, hjb.sig{j} = zeros(G.J, 1); end

end

