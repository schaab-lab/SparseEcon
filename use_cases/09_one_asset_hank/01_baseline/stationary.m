function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
r = x(1);
Y = x(2);

mc = (param.epsilon - 1) / param.epsilon;
w  = mc;
Pi = (1-mc) * Y;
N  = Y;

tau = Pi + param.tau_lab * w * N - param.G - r*param.gov_bond_supply;


%% VFI
G.r = r; G.w = w; G.tau = tau;

% Boundary conditions:
c_guess = r * G.a + (1-param.tau_lab) * w * param.zz + tau;

f = @(x,a) x - tau - r*a - ((1-param.tau_lab)*w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta);
f_prime = @(x) 1 + param.gamma/param.eta*((1-param.tau_lab)*w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta-1);

[G.c_right, ~] = newton_nonlin(f, f_prime, c_guess, param.amax, param.crit);
[G.c_left,  ~] = newton_nonlin(f, f_prime, c_guess, param.amin, param.crit);
[G.c0, ~]      = newton_nonlin(f, f_prime, c_guess,        G.a, param.crit);

% State-constrained boundary conditions:
left_bound  = param.u1(G.c_left);
right_bound = param.u1(G.c_right);

for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) sparse_project(left_bound(:, j),  points, G);
    BC{1}.right.f = @(points) sparse_project(right_bound(:, j), points, G);
    G = gen_FD(G, BC, num2str(j));
    G_dense = gen_FD(G_dense, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G,'V0'), G.V0 = param.u(c_guess) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, mu_dense{j}  = G.BH_dense * hjb.mu{j}; end
for j = 1:param.discrete_types, sig_dense{j} = G.BH_dense * hjb.sig{j}; end

g = KF(mu_dense, sig_dense, G_dense, param);
mass = sum(g .* G_dense.dx);

%% MARKET CLEARING
BH = sum(sum( G_dense.a .* g .* G_dense.dx));
NH = sum(sum( (G.BH_dense * (param.zz .* hjb.l)) .* g .* G_dense.dx));
C  = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S  = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

excess_bonds   = BH - param.gov_bond_supply;
excess_saving  = S;
excess_supply  = Y - C;
excess_labor   = N - NH;

diff = [excess_bonds, excess_supply];

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s; ss.mass = mass;
ss.N = N; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; ss.Pi = Pi; ss.tau = tau;
ss.excess_supply = excess_supply; ss.excess_bonds = excess_bonds; 
ss.excess_labor = excess_labor; ss.excess_saving = excess_saving;

end




