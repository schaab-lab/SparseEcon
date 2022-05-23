function [diff, G, G_dense, ss] = stationary(x, piw, G, G_dense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end
K = x(2);
N = x(3);

pi = piw;
i  = r + pi;

Z  = param.Z;
Y  = Z * K.^param.alpha .* N.^(1-param.alpha);
rk = param.alpha * Y / K;
w  = (1-param.alpha)  * Y / N;
I  = param.delta * K;
tau = param.tau_lab * w * N - param.G - r*param.gov_bond_supply;

ltau  = 15;
ltau0 = rk * (param.kmax*0.999)^(1-ltau);
xi    = param.xi * ltau0 * G.k .^ ltau;


%% VFI
G.income_a = r * G.a + rk * G.k - xi + (1-param.tau_lab) * w .* param.zz .* N + tau;
G.income_k = - param.delta * G.k;

G.N = N; G.piw = piw;

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
if ~isfield(G, 'V0'), G.V0 = (param.u(G.income_a) - param.v(G.N)) / param.rho; end

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
S = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

u1z = param.zz .* (G.BH_dense * param.u1(hjb.c));
M = sum(sum(u1z .* g .* G_dense.dx));

IH  = sum(sum( (G.BH_dense * hjb.iota) .* g .* G_dense.dx));
KH  = sum(sum( G_dense.k .* g .* G_dense.dx));
Chi = sum(sum( adjcostfn(G.BH_dense * hjb.iota, G_dense.k, param) .* g .* G_dense.dx));
Xi  = sum(sum( (G.BH_dense * xi) .* g .* G_dense.dx));

excess_bonds   = param.gov_bond_supply - B;
excess_saving  = S;
excess_capital = KH - K;
excess_invest  = IH - I;
excess_labor   = param.v1(N) - (param.epsilon - 1)/param.epsilon * (1+param.tau_L) * (1-param.tau_lab) * w * M;
excess_goods   = Y - C - I - Chi - Xi - param.G;

diff = [excess_bonds, excess_capital, excess_labor]';

ss.V = V; ss.g = g; ss.iota = hjb.iota; ss.w = w; ss.A = hjb.A; ss.m = hjb.m; ss.c = hjb.c; ss.s = hjb.s; ss.c0 = hjb.c0;
ss.mass = mass; ss.excess_goods = excess_goods; ss.excess_bonds = excess_bonds; ss.excess_capital = excess_capital; 
ss.excess_labor = excess_labor; ss.excess_invest = excess_invest; ss.excess_saving = excess_saving;
ss.Y = Y; ss.Z = Z; ss.C = C; ss.K = K; ss.KH = KH; ss.B = B; ss.I = I; ss.M = M; ss.u1z = u1z; ss.piw = piw; ss.pi = pi;
ss.r = r; ss.rk = rk; ss.S = S; ss.N = N; ss.Xi = Xi; ss.tau = tau; ss.Chi = Chi;
switch param.shock_type
    case 'TFP',       ss.shock = Z; 
    case 'demand',    ss.shock = param.rho; 
    case 'cost-push', ss.shock = param.epsilon; 
end

end




