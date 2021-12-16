function hjb = HJB(V, G, param)

VaF = deriv_sparse(G, V, 1, 'D1F');
VaB = deriv_sparse(G, V, 1, 'D1B');

cF = param.u1inv(VaF);
cB = param.u1inv(VaB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

IF = (sF>0) .* (G.grid(:,1)<1);  
IB = (sB<0) .* (IF==0) .* (G.grid(:,1)>0);
I0 = (1-IF-IB);

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

% COLLECT OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

% for j = 1:param.discrete_types, hjb.mu{j}  = s(:,j); end
% for j = 1:param.discrete_types, hjb.sig{j} = zeros(G.J, 1); end

end

