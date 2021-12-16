function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end

K = (param.alpha * exp(param.Zmean) / (r + param.delta)) ^ (1/(1-param.alpha)) * param.L;
w = (1-param.alpha) * exp(param.Zmean) * K.^param.alpha .* param.L.^(-param.alpha);
Y = exp(param.Zmean) * K^param.alpha * param.L^(1-param.alpha);


%% VFI
G.income = r * G.a + w .* param.zz;

% Boundary conditions:
left_bound  = param.u1(G.income(G.grid(:, 1) == 0, :));
right_bound = param.u1(G.income(G.grid(:, 1) == 1, :));
for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end

% Initialize guess V0
if ~isfield(G,'V0'), G.V0 = param.u(G.income) / param.rho; end

% Solve VFI
[V, hjb] = VFI_DSS(G, param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, mu_dense{j}  = G.BH_dense * hjb.mu{j}; end
for j = 1:param.discrete_types, sig_dense{j} = G.BH_dense * hjb.sig{j}; end

g = KF(mu_dense, sig_dense, G_dense, param);


%% MARKET CLEARING
KH = sum(sum( G_dense.a .* g * G_dense.dx));
C  = sum(sum( (G.BH_dense * hjb.c) .* g * G_dense.dx));
S  = sum(sum( (G.BH_dense * hjb.s) .* g * G_dense.dx));

excess_supply  = Y - C - param.delta*KH;
excess_capital = K - KH;

diff = excess_capital;

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s; ss.mass = sum(g * G_dense.dx);
ss.K = K; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; 
ss.excess_supply = excess_supply; ss.excess_capital = excess_capital; 

end


