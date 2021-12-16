function hjb = HJB(V, G, param)


%% VF DERIVATIVES
num0 = 1e-8; % numerical 0 for upwind scheme

VkF = deriv_sparse(G, V, 1, 'D1F');
VkB = deriv_sparse(G, V, 1, 'D1B');

VkF = max(VkF, num0);
VkB = max(VkB, num0);

%% UPWIND
cF = param.u1inv(VkF);
cB = param.u1inv(VkB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

HF = param.u(cF) + VkF .* sF;
HB = param.u(cB) + VkB .* sB;

IF = (sF > num0) .* ((sB > num0) + (sB <-num0).*(HF>=HB));
IB = (sB <-num0) .* ((sF <-num0) + (sF > num0).*(HB>=HF));
I0 = ~IF & ~IB;

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

%% OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

end

