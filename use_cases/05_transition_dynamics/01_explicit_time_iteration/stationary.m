function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
K = x(1);
L = x(2); % with inelastic labor supply: =param.L

Y  = K^param.alpha * L^(1-param.alpha);
rk = param.alpha * Y / K;
w  = (1-param.alpha) * Y / L;
r  = rk - param.delta;


%% VFI
G.r = r;
G.w = w;

% Boundary conditions:
c_guess = r * G.k + w * param.zz;

f = @(x,k) x - r*k - (w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta);
f_prime = @(x) 1 + param.gamma/param.eta*(w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta-1);

[G.c_right,~] = newton_nonlin(f, f_prime, c_guess, param.kmax, param.crit);
[G.c_left, ~] = newton_nonlin(f, f_prime, c_guess, param.kmin, param.crit);
[G.c0,~]      = newton_nonlin(f, f_prime, c_guess,        G.k, param.crit);

% State-constrained boundary conditions:
left_bound  = param.u1(G.c_left( G.grid(:, 1) == 0, :));
right_bound = param.u1(G.c_right(G.grid(:, 1) == 1, :));

for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j)  * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = param.u(c_guess) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, mu_dense{j}  = G.BH_dense * hjb.mu{j}; end
for j = 1:param.discrete_types, sig_dense{j} = G.BH_dense * hjb.sig{j}; end

g = KF(mu_dense, sig_dense, G_dense, param);


%% MARKET CLEARING
KH = sum(sum( G_dense.k .* g .* G_dense.dx));
LH = sum(sum( (G.BH_dense * (param.zz .* hjb.l)) .* g .* G_dense.dx));
C  = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S  = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

excess_supply  = Y - C - param.delta*KH;
excess_capital = K - KH;
excess_labor   = L - LH;

diff = [excess_capital, excess_labor];

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.K = K; ss.L = L; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; 
ss.excess_supply = excess_supply; ss.excess_capital = excess_capital; ss.excess_labor = excess_labor;

end




