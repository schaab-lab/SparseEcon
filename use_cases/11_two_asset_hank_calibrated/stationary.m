function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end
K = x(2); 
L = x(3);


%% AGGREGATES
mc = (param.epsilonF - 1) / param.epsilonF;

Z  = param.Z;
Y  = Z * K.^param.alpha .* L.^(1-param.alpha);
U  = (1-param.L);
H  = L / (1-U);
rk = mc * param.alpha * Y / K;
w  = mc * (1-param.alpha)  * Y / L;
Pi = (1-mc) * Y;

tau = Pi + param.tau_lab * w * L - param.UI * U - param.G - r*param.gov_bond_supply;

gross_total_capital_accumulation = param.delta * K;
Q = param.solve_for_Q_from_cap_accumulation(gross_total_capital_accumulation, K);

gross_total_investment_expenditure = param.gross_total_investment_expenditure(Q, K, 0);
I   = gross_total_investment_expenditure;
PiQ = param.PiQ(Q, K, param.ZQmean);
muK = param.gross_total_capital_accumulation(Q, K, param.ZQmean) - param.delta * K;

cap_income_liquid   = rk + PiQ / K;
cap_income_illiquid = param.deathrate - param.delta;

ltau  = 15;
ltau0 = cap_income_liquid * (param.kmax*0.999)^(1-ltau);
xi    = param.xi * ltau0 * G.k .^ ltau;


%% VFI
G.income_a = (r + param.deathrate) .* G.a ...
             + cap_income_liquid .* G.k - xi ...
             + (1-param.tau_lab) * w .* param.zz .* H + tau + param.UI .* [1, 0];
G.income_k = cap_income_illiquid .* G.k;

G.Q = Q; G.H = H;

% State-constrained boundary conditions:
% left_bound  = param.u1(G.income_a);
% right_bound = param.u1(G.income_a);
% 
% for j = 1:param.discrete_types
%     BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
%     BC{1}.left.f  = @(points) sparseProject(left_bound(:, j),  points, G);
%     BC{1}.right.f = @(points) sparseProject(right_bound(:, j), points, G);
%     BC{2}.left.type = '0'; BC{2}.right.type = '0';
%     G = gen_FD(G, BC, num2str(j));
%     GDense = gen_FD(GDense, BC, num2str(j));
% end

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = (param.u(G.income_a) - param.v(G.H)) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, param);
if isnan(V), diff = NaN(1); return; end


%% OUTPUT VF as next guess
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, sc_dense{j} = G.BH_dense * hjb.sc(:, j);   end
for j = 1:param.discrete_types, si_dense{j} = G.BH_dense * hjb.si(:, j);   end
for j = 1:param.discrete_types, mi_dense{j} = G.BH_dense * hjb.iota(:, j); end
for j = 1:param.discrete_types, mk_dense{j} = G.BH_dense * G.income_k;    end

g = KF(sc_dense, si_dense, mi_dense, mk_dense, G_dense, param);

mass = sum(g * G_dense.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end


%% MARKET CLEARING
B = sum(sum( G_dense.a .* g .* G_dense.dx));
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = - param.deathrate * B + sum(sum( (G.BH_dense * hjb.s).*g .* G_dense.dx));

IH  = sum(sum( (G.BH_dense * hjb.iota) .* g .* G_dense.dx));
KH  = sum(sum( G_dense.k .* g .* G_dense.dx));
Chi = sum(sum( adjcostfn(G.BH_dense * hjb.iota, G_dense.k, param) .* g .* G_dense.dx));
Xi  = sum(sum( (G.BH_dense * xi) .* g .* G_dense.dx));

M1 = muK;
M2 = - param.deathrate * KH + sum(sum( (G.BH_dense * hjb.m) .*g  .* G_dense.dx));

Lambda = sum(sum( ( param.zz .* param.u1(G.BH_dense * hjb.c) ) .* g .* G_dense.dx));

excess_bonds   = param.gov_bond_supply - B;
excess_saving  = S;
excess_capital = KH - K;
excess_labor   = param.v1(H) - (param.epsilonW - 1)/param.epsilonW * (1-param.tau_lab) * w * Lambda;

excess_cap_production = param.gross_total_capital_accumulation(Q, K, param.ZQmean) - IH;
excess_goods = Y - C - param.gross_total_investment_expenditure(Q, K, param.ZQmean) - Chi - Xi - param.G;

zBar = sum(sum( param.zz .* g .* G_dense.dx));

diff = [excess_bonds, excess_capital, excess_labor]';

ss.V = V; ss.g = g; ss.iota = hjb.iota; ss.w = w; ss.A = hjb.A; ss.m = hjb.m; ss.c = hjb.c; ss.s = hjb.s; ss.c0 = hjb.c0;
ss.mass = mass; ss.excess_goods = excess_goods; ss.excess_bonds = excess_bonds; ss.excess_capital = excess_capital; 
ss.excess_labor = excess_labor; ss.excess_cap_production = excess_cap_production; ss.excess_saving = excess_saving;
ss.Y = Y; ss.C = C; ss.K = K; ss.KH = KH; ss.B = B; ss.Xi = Xi; ss.I = I;
ss.r = r; ss.rk = rk; ss.S = S; ss.L = L; ss.H = H; ss.zBar = zBar; ss.Pi = Pi; ss.Xi=Xi;
ss.M1 = M1; ss.M2 = M2; ss.tau = tau; ss.Q = Q; ss.Chi = Chi; ss.Pi = Pi; ss.PiQ = PiQ;
switch param.shock_type
    case 'TFP',       ss.shock = Z; 
    case 'demand',    ss.shock = param.rho; 
    case 'cost-push', ss.shock = param.epsilon; 
end

end




