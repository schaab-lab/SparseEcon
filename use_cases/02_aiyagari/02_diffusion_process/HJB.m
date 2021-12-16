function hjb = HJB(V, G, param)

num0 = 1e-8; % numerical 0 for upwind scheme

VkF = deriv_sparse(G, V, 1, 'D1F');
VkB = deriv_sparse(G, V, 1, 'D1B');

cF = param.u1inv(VkF);
cB = param.u1inv(VkB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
I0 = ~IF & ~IB;

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

% COLLECT OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

end

