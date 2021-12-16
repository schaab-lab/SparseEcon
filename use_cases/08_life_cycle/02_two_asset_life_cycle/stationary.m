function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

K = x(1); 
L = x(2);
Q = x(3);
tau = x(4);

%% AGGREGATES
r  = param.r; 
Y  = K.^param.alpha .* L.^(1-param.alpha);
H  = L / param.L;
rk = param.alpha * Y / K;
w  = (1-param.alpha)  * Y / L;

PiQ = param.PiQ(Q, K, param.ZQmean);
muK = param.gross_total_capital_accumulation(Q, K, param.ZQmean) - param.delta * K;

cap_income_liquid   = rk + PiQ / K;
cap_income_illiquid = param.deathrate - param.delta;

ltau  = 15;
ltau0 = cap_income_liquid * (param.kmax*0.999)^(1-ltau);
xi    = param.xi * ltau0 * G.k .^ ltau;


%% VFI
G.income_a = (r + param.deathrate) .* G.a + cap_income_liquid .* G.k - xi ...
             + (1-param.tau_lab) * w .* param.zz .* H + tau + param.UI .* [1, 0];
G.income_k = cap_income_illiquid .* G.k;

G.Q = Q; G.H = H;

% State-constrained boundary conditions:
% left_bound  = param.u1(G.income_a);
% right_bound = param.u1(G.income_a);

right_bound_T = 1e-8 * param.u(1e-5 + G.a) .* [1,1];
% right_bound_T(G.a<0,:) = - 100 * G.a(G.a<0)/param.amin .* [1,1];

for j = 1:param.discrete_types
    %BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    %BC{1}.left.f  = @(points) sparseProject(left_bound(:,j),  points, G);
    %BC{1}.right.f = @(points) sparseProject(right_bound(:,j), points, G);
    BC{1}.left.type = '0'; BC{1}.right.type = '0';
    BC{2}.left.type = '0'; BC{2}.right.type = '0';
    BC{3}.left.type = '0'; BC{3}.right.type = 'D';
    BC{3}.right.f = @(points) sparse_project(right_bound_T(:,j), points, G);
    G = gen_FD_dirichlet(G, BC, num2str(j));
    BC{3}.right.type = 'D_birth';
    G_dense = gen_FD_dirichlet(G_dense, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G,'V0')
    G.V0 = param.u(G.income_a);
    G.V0(G.t==param.tmax, :) = right_bound_T(G.t==param.tmax, :);
end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF as next guess
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, sc_dense{j} = G.BH_dense * hjb.sc(:, j);   end
for j = 1:param.discrete_types, si_dense{j} = G.BH_dense * hjb.si(:, j);   end
for j = 1:param.discrete_types, mi_dense{j} = G.BH_dense * hjb.iota(:, j); end
for j = 1:param.discrete_types, mk_dense{j} = G.BH_dense * G.income_k;     end

g = KF(sc_dense, si_dense, mi_dense, mk_dense, G_dense, param);

mass = sum(g * G_dense.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end


%% MARKET CLEARING
B = sum(sum( G_dense.a .* g .* G_dense.dx));
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = - param.deathrate * B + sum(sum( (G.BH_dense * hjb.s).*g .* G_dense.dx));

IH  = sum(sum( (G.BH_dense * hjb.iota) .* g .* G_dense.dx));
KH  = sum(sum( G_dense.k .* g .* G_dense.dx));
Chi = sum(sum( adjcostfn(G.BH_dense * hjb.iota,G_dense.k,param) .* g .* G_dense.dx));
Xi  = sum(sum( (G.BH_dense * xi) .* g .* G_dense.dx));

M1 = muK;
M2 = - param.deathrate * KH + sum(sum( (G.BH_dense * hjb.m) .*g  .* G_dense.dx));

Lambda = sum(sum( ( param.zz .* param.u1(G.BH_dense * hjb.c) ) .* g .* G_dense.dx));

I = param.gross_total_investment_expenditure(Q, K, param.ZQmean);

excess_bonds   = param.gov_bond_supply - B;
excess_capital = KH - K;
excess_labor   = param.v1(H) - (param.epsilonW - 1)/param.epsilonW * (1-param.tau_lab) * w * Lambda;
excess_tax     = tau - (param.tau_lab * w * L - param.UI * mass(1) - param.G - r*param.gov_bond_supply);

excess_cap_production = param.gross_total_capital_accumulation(Q, K, param.ZQmean) - IH;
excess_goods = Y - C - param.gross_total_investment_expenditure(Q, K, param.ZQmean) - Chi - Xi - param.G;

zBar = sum(sum( param.zz .* g .* G_dense.dx));

diff = [excess_capital, excess_labor, excess_cap_production, excess_tax]';

ss.V = V; ss.g = g; ss.iota = hjb.iota; ss.w = w; ss.A = hjb.A; ss.m = hjb.m; ss.c = hjb.c; ss.s = hjb.s; ss.c0 = hjb.c0;
ss.mass = mass; ss.excess_goods = excess_goods; ss.excess_bonds = excess_bonds; ss.excess_capital = excess_capital; 
ss.excess_labor = excess_labor; ss.excess_cap_production = excess_cap_production;
ss.Y = Y; ss.C = C; ss.K = K; ss.KH = KH; ss.B = B; ss.Xi = Xi; ss.I = I;
ss.r = r; ss.rk = rk; ss.S = S; ss.L = L; ss.H = H; ss.zBar = zBar; ss.Xi=Xi;
ss.M1 = M1; ss.M2 = M2; ss.tau = tau; ss.Q = Q; ss.Chi = Chi; ss.PiQ = PiQ;

end




