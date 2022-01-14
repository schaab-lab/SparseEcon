function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end

Y = exp(param.shock_mean) * param.L;
w = exp(param.shock_mean);


%% VFI
G.income = r * G.a + w .* G.z;

% State-constrained boundary conditions:
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);

BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
BC{1}.left.f  = @(points) sparse_project(left_bound,  points, G);
BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
BC{2}.left.type = '0'; BC{2}.right.type = '0';
G = gen_FD(G, BC);
G_dense = gen_FD(G_dense, BC); % Note this is not actually necessary for Huggett

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = param.u(G.income) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
s_dense = G.BH_dense * hjb.s;
g = KF(s_dense, G_dense, param);


%% MARKET CLEARING
B = sum(sum( G_dense.a .* g .* G_dense.dx));
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));
L = sum(sum( G_dense.z .* g .* G_dense.dx));

if abs(L - param.L) > 1e-8, error('Aggregate labor not normalized.\n'); end

diff = B;

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.B = B; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; ss.excess_supply = Y - C;


end




