function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
K   = x(1);
tau = x(2);

Y  = K^param.alpha * param.L^(1-param.alpha);
rk = param.alpha * Y / K;
w  = (1-param.alpha) * Y / param.L;
r  = rk - param.delta;
rb = r + param.wedge;

%% VFI

G.income = (r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + w * G.z + tau;

% State-constrained boundary conditions:
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);

BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
BC{1}.left.f  = @(points) sparse_project(left_bound , points, G);
BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
BC{2}.left.type = '0'; BC{2}.right.type = '0';
G = gen_FD(G, BC);
G_dense = gen_FD(G_dense, BC);

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = param.u(G.income) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
s_dense = G.BH_dense * hjb.s;
g = KF(s_dense, G_dense, param);
if any(isnan(g)), diff = NaN(1); return; end


%% MARKET CLEARING
AH = sum(sum( G_dense.a .* g .* G_dense.dx));
C  = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S  = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));
KH = AH - param.gov_bond_supply;

excess_capital = K - KH;
excess_supply  = Y - C - param.delta*KH;
excess_rebate  = tau - param.wedge * sum(sum( G_dense.a .* g .* (G_dense.a<0) .* G_dense.dx));
% excess_rebate  = tau - param.wedge * sum( ((G.BHDense * (G.a.*(G.a<0))) .*g))*GDense.dx;

diff = [excess_capital, excess_rebate];

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.K = K; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; ss.tau = tau;
ss.excess_supply = excess_supply; ss.excess_capital = excess_capital;

%{
max(abs( ss.s - (G.income - ss.c) ))
max(abs( ss.s - ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + w * G.z + tau - ss.c) ))

max(abs( sum((G.BHDense * hjb.s).*g)*GDense.dx ...
    - sum( (G.BHDense * ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + w * G.z + tau - ss.c)).*g)*GDense.dx ))

S - sum( (G.BHDense * ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + w * G.z + tau - ss.c)).*g)*GDense.dx

S + C - sum( (G.BHDense * ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + w * G.z + tau)).*g)*GDense.dx
S + C - w*param.L - sum( (G.BHDense * ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0) + tau)).*g)*GDense.dx
S + C - w*param.L - tau - sum( (G.BHDense * ((r * G.a).*(G.a>=0) + (rb * G.a).*(G.a<0))).*g)*GDense.dx
S + C - w*param.L - tau - ...
    sum( (G.BHDense * (((rk - param.delta) * G.a).*(G.a>=0) + ((rk - param.delta+param.wedge) * G.a).*(G.a<0))).*g)*GDense.dx

S + C - w*param.L - tau - (rk-param.delta)*AH - param.wedge * sum( ((G.BHDense * (G.a.*(G.a<0))) .*g))*GDense.dx

S + C - w*param.L - tau - (rk-param.delta)*AH - param.wedge * sum( GDense.a.*(GDense.a<0) .*g)*GDense.dx

%}

end




