function hjb = HJB(V, G, param)

num0 = 1e-8; % numerical 0 for upwind scheme

VaF = zeros(G.J, param.discrete_types);
VaB = zeros(G.J, param.discrete_types);
for j = 1:param.discrete_types
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

cF = param.u1inv(VaF);
cB = param.u1inv(VaB);

lF = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(cF));
lB = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(cB));
l0 = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(G.c0));

sF = G.r.*G.a + (1-param.tau_lab)*G.w.*param.zz.*lF + G.tau - cF;
sB = G.r.*G.a + (1-param.tau_lab)*G.w.*param.zz.*lB + G.tau - cB;
s0 = G.r.*G.a + (1-param.tau_lab)*G.w.*param.zz.*l0 + G.tau - G.c0;

IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
I0 = ~IF & ~IB;

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + G.c0.*I0;
l = lF.*IF + lB.*IB + l0.*I0;
u = param.u(c) - param.v(l);

% COLLECT OUTPUT
hjb.c = c; hjb.l = l; hjb.s = s; hjb.u = u; hjb.sF = sF; hjb.sB = sB;

for j = 1:param.discrete_types, hjb.mu{j}  = s(:, j); end
for j = 1:param.discrete_types, hjb.sig{j} = zeros(G.J, 1); end

% TESTS
assert(all(all( abs(s0) < 1e-10 )));
assert(all(all( s(G.grid(:, 1) == 1,:) <= 0 )));
assert(all(all( s(G.grid(:, 1) == 0,:) >= 0 )));
assert(all(all( abs(l - ((1-param.tau_lab)*G.w * param.zz .* c.^(-param.gamma)).^(1/param.eta)) < 1e-10 )));

assert(all(all( abs(cF(G.grid(:, 1) == 1, :) - G.c_right) < 1e-10 )));
assert(all(all( abs(cB(G.grid(:, 1) == 0, :) - G.c_left ) < 1e-10 )));

end

