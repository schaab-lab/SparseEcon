function hjb = HJB(V, G, param)


%% VF DERIVATIVES
num0 = 1e-8; % numerical 0 for upwind scheme

[VaF, VaB, VkF, VkB] = deal(zeros(G.J, param.discrete_types));
for j = 1:param.discrete_types
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    VkF(:, j) = deriv_sparse(G, V(:, j), 2, 'D1F', num2str(j));
    VkB(:, j) = deriv_sparse(G, V(:, j), 2, 'D1B', num2str(j));
end

VkF = max(VkF, num0);
VkB = max(VkB, num0);
VaF = max(VaF, num0);
VaB = max(VaB, num0);


%% CONSUMPTION + SAVINGS
cF = param.u1inv(VaF); cF(G.a == param.amax, :) = G.income_a(G.a == param.amax, :);
cB = param.u1inv(VaB); cB(G.a == param.amin, :) = G.income_a(G.a == param.amin, :);
c0 = G.income_a;

scF = G.income_a - cF;
scB = G.income_a - cB;

IF = (scF > num0);
IB = (scB <-num0) & ~IF; 
I0 = ~IF & ~IB;

c = cF.*IF + cB.*IB + c0.*I0;
sc = scF.*IF + scB.*IB;


%% INVESTMENT
iotaFF = adjcostfn1inv(VkF./VaF-1, G.k, param); iotaFF(G.k == param.kmax, :) = 0;
iotaFB = adjcostfn1inv(VkB./VaF-1, G.k, param); iotaFB(G.k == param.kmin, :) = 0;
iotaBF = adjcostfn1inv(VkF./VaB-1, G.k, param); iotaBF(G.k == param.kmax, :) = 0;
iotaBB = adjcostfn1inv(VkB./VaB-1, G.k, param); iotaBB(G.k == param.kmin, :) = 0;

siFF = - iotaFF - adjcostfn(iotaFF, G.k, param); siFF(G.a == param.amax, :) = 0;
siFB = - iotaFB - adjcostfn(iotaFB, G.k, param); siFB(G.a == param.amax, :) = 0;
siBF = - iotaBF - adjcostfn(iotaBF, G.k, param); siBF(G.a == param.amin, :) = 0;
siBB = - iotaBB - adjcostfn(iotaBB, G.k, param); siBB(G.a == param.amin, :) = 0;

IFF = (siFF > num0) & (iotaFF > num0);
IFB = (siFB > num0) & (iotaFB <-num0) & ~IFF;
IBF = (siBF <-num0) & (iotaBF > num0) & ~IFF & ~IFB;
IBB = (siBB <-num0) & (iotaBB <-num0) & ~IFF & ~IFB & ~IBF;

iota = iotaFF.*IFF + iotaFB.*IFB + iotaBF.*IBF + iotaBB.*IBB;
si = - iota - adjcostfn(iota, G.k, param);


%% OUTPUT
u = param.u(c) - param.v(G.N) - param.delta/2 * G.piw^2;
s = sc + si;
m = G.income_k + iota;

hjb.c = c; hjb.s = s; hjb.m = m; hjb.iota = iota; hjb.sc = sc; hjb.si = si; hjb.u = u; hjb.c0 = c0;


end


