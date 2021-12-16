function hjb = HJB(V, G, param)

num0 = 1e-10;

VaF = zeros(G.J, param.discrete_types);
VaB = zeros(G.J, param.discrete_types);
for j = 1:param.discrete_types
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

VaF = max(VaF, num0);
VaB = max(VaB, num0);

cF = param.u1inv(VaF);
cB = param.u1inv(VaB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

IF = (sF > num0);
IB = (sB <-num0) .* ~IF;
I0 = (1-IF-IB);

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

% COLLECT OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

% for j = 1:param.discrete_types, hjb.mu{j}  = s(:,j); end
% for j = 1:param.discrete_types, hjb.sig{j} = zeros(G.J, 1); end

end

