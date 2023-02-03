function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end
N = x(2);

Z = param.Z;
Y = Z * N;
w = Z;


%% VFI
G.income = r * G.a + w .* param.zz .* N;
G.N = N;

% State-constrained boundary conditions:
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);

for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) sparse_project(left_bound(:, j),  points, G);
    BC{1}.right.f = @(points) sparse_project(right_bound(:, j), points, G);
    BC{2}.left.type = '0'; BC{2}.right.type = '0';
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = (param.u(G.income) - param.v(G.N)) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, param); 
if isnan(V), diff = NaN(1); return; end


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
s_dense = G.BH_dense * hjb.s;
g = KF(s_dense, G_dense, param);


%% MARKET CLEARING
B = sum(sum( G_dense.a .* g .* G_dense.dx));
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

m = param.zz .* (G.BH_dense * param.u1(hjb.c));
M = sum(sum(m .* g .* G_dense.dx));

ss.excess_supply = Y - C;
ss.excess_saving = S;
ss.excess_bonds = B;
ss.excess_union = (param.epsilon-1)/param.epsilon * (1+param.tau_L) * param.Z * M - param.v1(N);

diff = [B, ss.excess_union];

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s; ss.m = m; ss.N = N; ss.M = M; ss.Z = Z;
ss.B = B; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w;
ss.A = hjb.A; ss.Aa = hjb.Aa; ss.Az = hjb.Az; ss.AT = hjb.A'; ss.u = hjb.u;
switch param.shock_type
    case 'TFP',       ss.shock = Z; 
    case 'demand',    ss.shock = param.rho; 
    case 'cost-push', ss.shock = param.epsilon; 
end

end




