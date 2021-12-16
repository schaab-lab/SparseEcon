function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end
K = x(2);


%% AGGREGATES
Y = K.^param.alpha .* param.L.^(1-param.alpha);

w = (1-param.alpha) .* Y ./ param.L;
D = param.alpha .* Y ./ K;

muR = D - param.delta;

if muR-r < 0.001, diff = NaN(1); return; end


%% VFI
G.r = r; G.w = w; G.muR = muR;
G.income = G.r .* G.n + G.w .* param.zz;

% State-constrained boundary conditions:
left_bound = param.u1(G.income);
for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'fxF=coef*fxx';
    BC{1}.left.f  = @(points) sparse_project(left_bound(:, j), points, G);
    BC{1}.right.f  = @(points) -param.nmax / param.gamma * ones(size(points, 1));
    G = gen_FD_chain_rule(G, BC, num2str(j)); % GDense not necessary because "inward-pointing"
end

% Initialize guess V0:
if ~isfield(G,'V0')
    % G.V0 = param.u(G.income + (G.muR-G.r).^2 .* G.n / (param.gamma*param.sigR2)) / param.rho; 
    G.V0 = param.u(G.income) / param.rho; 
end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, mu_dense{j}  = G.BH_dense * hjb.s(:, j); end
for j = 1:param.discrete_types, sig_dense{j} = G.BH_dense * (param.sigR * hjb.theta(:, j) .* G.n); end

g = KF(mu_dense, sig_dense, G_dense, param);


%% MARKET CLEARING
N = sum(sum( G_dense.n .* g .* G_dense.dx));
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

B = sum(sum( (1-G.BH_dense * hjb.theta) .* G_dense.n .* g .* G_dense.dx));
KS = sum(sum( (G.BH_dense * hjb.theta) .* G_dense.n .* g .* G_dense.dx));

excess_supply  = Y - C - param.delta * K;
excess_capital = K - KS;
excess_savings = S;
excess_bonds   = B;

diff = [B, excess_supply]';

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s; ss.theta = hjb.theta;
ss.B = N; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; ss.K = K;
ss.excess_capital = excess_capital; ss.excess_supply = excess_supply; 
ss.excess_savings = excess_savings; ss.excess_bonds = excess_bonds;

end




