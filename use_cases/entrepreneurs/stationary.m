function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end

% Public firm capital to labor ratio xi = K/L & wage
xi = (param.eta*param.Aprod*param.Bprod/(r+param.delta)).^(1/(1-param.eta));
w = (1-param.eta)*param.Aprod*param.Bprod*xi.^param.eta;
 
%% FIXED COSTS IN UNITS OF CAPITAL
% Optimal capital and labor choices corresponding to the unproductive technology
kuU = (G.z.*param.Aprod.*param.BU).^(1/(1-param.alpha-param.beta)).*(param.alpha/(r+param.delta)).^((1-param.beta)/(1-param.alpha-param.beta)).*(param.beta/w).^(param.beta/(1-param.alpha-param.beta)) + param.fkU;
kU  = max(min(param.lambda*G.a, kuU), 0);
lU  = (param.beta*G.z*param.Aprod*param.BU./w).^(1/(1-param.beta)).*kU.^(param.alpha/(1-param.beta)) + param.flU;

% Optimal capital and labor choices corresponding to the productive technology
kuP = (G.z.*param.Aprod.*param.BP).^(1/(1-param.alpha-param.beta)).*(param.alpha/(r+param.delta)).^((1-param.beta)/(1-param.alpha-param.beta)).*(param.beta/w).^(param.beta/(1-param.alpha-param.beta)) + param.fkP;
kP  = max(min(param.lambda*G.a, kuP), 0);
lP  = (param.beta*G.z*param.Aprod*param.BP./w).^(1/(1-param.beta)).*kP.^(param.alpha/(1-param.beta)) + param.flP;

% Income maximization: occupation and technology choice
YU  = G.z*param.Aprod*param.BU.*max(kU-param.fkU,0).^param.alpha.*max(lU-param.flU, 0).^param.beta;
YP  = G.z*param.Aprod*param.BP.*max(kP-param.fkP,0).^param.alpha.*max(lP-param.flP, 0).^param.beta;
PiU = YU - (r+param.delta).*kU - w.*lU - param.fyU;
PiP = YP - (r+param.delta).*kP - w.*lP - param.fyP;
M  = max(w.*G.z.^param.theta, max(PiU,PiP));

% Indicators of the 3 different occupation and technology choices
worker       = (w.*G.z.^param.theta > max(PiU, PiP)); 
unproductive = (PiU > max(w.*G.z.^param.theta, PiP)); 
productive   = (PiP >= max(w.*G.z.^param.theta, PiU)); 


%% VFI
G.income = M + r.*G.a;

% State-constrained boundary conditions:
left_bound  = param.u1(w*G.z.^param.theta + r*param.amin);
right_bound = param.u1(w*G.z.^param.theta + r*param.amax);

BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
BC{1}.left.f  = @(points) sparse_project(left_bound,  points, G);
BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
BC{2}.left.type = '0'; BC{2}.right.type = '0';
G = gen_FD(G, BC);
G_dense = gen_FD(G_dense, BC); % Note this is not actually necessary for Huggett

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = param.u(w*G.z + r*G.a) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, [], param);

%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
s_dense = G.BH_dense * hjb.s;
g = KF(s_dense, G_dense, param);

%% MARKET CLEARING
% capital
K_s = sum(sum( G_dense.a .* g .* G_dense.dx));
K_d = sum(sum( (G.BH_dense * (kU .* unproductive)) .* g .* G_dense.dx)) + sum(sum( (G.BH_dense * (kP .* productive)) .* g .* G_dense.dx));

% labor
L_s = sum(sum( (G.BH_dense * (G.z.^param.theta.*worker)) .*g .*G_dense.dx));
L_d = sum(sum( (G.BH_dense * (lU .* unproductive)) .* g .*G_dense.dx)) + sum(sum( (G.BH_dense * (lP .* productive)) .* g .*G_dense.dx));

% fraction of worker & (prod / unprod - entrepreneur)
frac_worker = sum(sum( (G.BH_dense * worker) .*g .*G_dense.dx));
frac_prod   = sum(sum( (G.BH_dense * productive) .*g .*G_dense.dx));
frac_unprod = sum(sum( (G.BH_dense * unproductive) .*g .*G_dense.dx));
frac_entrepreneur = frac_prod + frac_unprod;

% aggregates
C = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

% public firm
K_c = max(K_s - K_d,1e-7); % i.e. imposing capital market clearing
L_c = K_c / xi;
Y_c = param.Aprod*param.Bprod*K_c^param.eta*L_c^(1-param.eta);

% entrepreneurs output
Y_U = sum(sum( (G.BH_dense * (YU .* unproductive)) .* g .*G_dense.dx));
Y_P = sum(sum( (G.BH_dense * (YP .* productive)) .* g .*G_dense.dx));

excess_capital = K_d + K_c - K_s;
excess_labor   = L_d + L_c - L_s;
excess_supply  = Y_c + Y_U + Y_P - C - param.delta*K_s;

diff = excess_labor;

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.K_s = K_s; ss.K_d = K_d; ss.K_c = K_c;
ss.L_s = L_s; ss.L_d = L_d; ss.L_c = L_c;
ss.Y_U = Y_U; ss.Y_P = Y_P; ss.Y_c = Y_c;
ss.C = C; ss.S = S; ss.r = r; ss.w = w;
ss.excess_labor = excess_labor; ss.excess_capital = excess_capital; ss.excess_supply = excess_supply;
ss.frac_worker = frac_worker; ss.frac_prod = frac_prod; ss.frac_unprod = frac_unprod; ss.frac_entrepreneur = frac_entrepreneur;
end




