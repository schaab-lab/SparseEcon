Q = sim.Q(1);
I = sim.I(1);
IH = sim.IH(1);
K = sim.K(1);
M = sim.M(1);

IQ = param.solveForIota(Q) * K
iotaQ = IQ / K;

assert(abs( I - (IQ + param.Phi(iotaQ)*K) ) < 1e-10)

assert(abs( Q - (1 + param.PhiPrime(iotaQ)) ) < 1e-10)

IH == param.grossTotalCapitalAccumulation(Q,K,0)

sim.K(2) == sim.K(1) + param.dt(1) * (sim.grossTotalCapitalAccumulation(1) - param.delta * sim.K(1))

assert(abs( M - ( IH - param.delta * sim.K(1) )) < 1e-10)

(sim.KH(2) - sim.KH(1))/ param.dt(1) + param.delta * sim.KH(1) = IH ??

sum(sum( GDense.k .* sim.g{2} * GDense.dx )) - sum(sum( GDense.k .* sim.g{1} * GDense.dx ))

(sim.KH(2) - sim.KH(1))/ param.dt(1) == sim.M(1) ??






dgdt = AT * [sim.g{n}(:,1); sim.g{n}(:,2)] + ...
         param.deathrate * ([birthID(:,1); birthID(:,2)] / GDense.dx - [sim.g{n}(:,1); sim.g{n}(:,2)]);
     
sum(sum( [GDense.k; GDense.k] .* dgdt )) * GDense.dx


%% TRANSITION DEBUGGING
%{
DEBUGGING IQ = IH:

max(abs( sim.Q - sim.Q2 ))
max(abs( sim.Q - sim.Q3 ))

max(abs( sim.I - (sim.IQ + sim.Phi) ))

max(abs( sim.IQ - sim.IH ))

max(abs( sim.IQ - (sim.dK + param.delta * sim.K) ))
max(abs( sim.IQ(1:end-1) - (diff(sim.K)./sim.dt(1:end-1) + param.delta * sim.K(1:end-1)) ))
max(abs( sim.IH(1:end-1) - (diff(sim.KH)./sim.dt(1:end-1)+ param.delta * sim.KH(1:end-1)) ))

max(abs( sim.M - sim.IH + param.delta*sim.KH ))

=> IH and M are simply not representing the drift of KH! This is about the KF!

%}


%{
s = sim.s{1};
c = sim.c{1};
m = sim.m{1};
iota = sim.iota{1};

sDense = G.BHDense * s;

n = 1;
capIncomeLiquid   = sim.rk(n) + sim.PiQ(n) / sim.K(n);
capIncomeIlliquid = param.deathrate - param.delta;

ltau  = 15;
ltau0 = capIncomeLiquid * (param.kmax*0.999)^(1-ltau);
sim.xi{n} = param.xi * ltau0 * G.k .^ ltau;

G.income_a = (sim.r(n) + param.deathrate) .* G.a ...
             + capIncomeLiquid .* G.k - sim.xi{n} ...
             + (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) + sim.tau(n) + param.UI .* [1, 0];
G.income_k = capIncomeIlliquid .* G.k;


sum(sum(abs( s - (G.income_a - c - sim.Q(1).*iota - adjcostfn(iota, G.k, param)) )))
sum(sum(abs( m - (G.income_k + iota) )))

sum(sum(abs( s - (...
             (sim.r(n) + param.deathrate) .* G.a ...
             + capIncomeLiquid .* G.k - sim.xi{n} ...
             + (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) + sim.tau(n) + param.UI .* [1, 0] ...
            - c - sim.Q(1).*iota - adjcostfn(iota, G.k, param)) )))

sum(sum(abs( G.BHDense*(s - (...
             (sim.r(n) + param.deathrate) .* G.a ...
             + capIncomeLiquid .* G.k - sim.xi{n} ...
             + (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) + sim.tau(n) + param.UI .* [1, 0] ...
            - c - sim.Q(1).*iota - adjcostfn(iota, G.k, param))) .* sim.g{1} .* GDense.dx )))

sum(sum(abs( G.BHDense*(s ...
             - (sim.r(n) + param.deathrate) .* G.a ...
             - capIncomeLiquid .* G.k + sim.xi{n} ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - sim.tau(n) - param.UI .* [1, 0] ...
             + c + sim.Q(1).*iota + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx )))

sum(sum( G.BHDense*(s - (sim.r(n) + param.deathrate) .* G.a ...
             - capIncomeLiquid .* G.k + sim.xi{n} ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - sim.tau(n) - param.UI .* [1, 0] ...
             + c + sim.Q(1).*iota + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))

sim.S(1) + param.deathrate*sim.B(1) + sum(sum( G.BHDense*( - (sim.r(n) + param.deathrate) .* G.a ...
             - capIncomeLiquid .* G.k + sim.xi{n} ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - sim.tau(n) - param.UI .* [1, 0] ...
             + c + sim.Q(1).*iota + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))

sim.S(1) + param.deathrate*sim.B(1) - sim.r(1)*sim.B(1) - param.deathrate * sim.B(1) + sum(sum( G.BHDense*(  ...
             - capIncomeLiquid .* G.k + sim.xi{n} ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - sim.tau(n) - param.UI .* [1, 0] ...
             + c + sim.Q(1).*iota + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))

sim.S(1) - sim.r(1)*sim.B(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) ...
             + sum(sum( G.BHDense*(  ...
             - capIncomeLiquid .* G.k ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - sim.tau(n) - param.UI .* [1, 0] ...
             + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))

sim.S(1) - sim.r(1)*sim.B(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) - sim.tau(1) ...
             + sum(sum( G.BHDense*(  ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) - param.UI .* [1, 0] ...
             + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))

sim.S(1) - sim.r(1)*sim.B(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) - sim.tau(1) ...
- param.UI * sim.U(1) + ...
             + sum(sum( G.BHDense*(  ...
             - (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) ...
             + adjcostfn(iota, G.k, param)) .* sim.g{1} .* GDense.dx ))


sim.S(1) - sim.r(1)*sim.B(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) - sim.tau(1) ...
- param.UI * sim.U(1) - (1-param.tau_lab)*sim.w(1)*sim.L(1) + sim.Chi(1)

sim.S(1) - sim.r(1)*sim.B(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) ...
- (sim.Pi(1) + param.tau_lab * sim.w(1) .* sim.L(1) - param.UI * sim.U(1) - sim.G(1) - sim.r(1) * param.gov_bond_supply) ...
- param.UI * sim.U(1) - (1-param.tau_lab)*sim.w(1)*sim.L(1) + sim.Chi(1)


sim.S(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) ...
- sim.Pi(1) - sim.w(1)*sim.L(1) + sim.Chi(1)

sim.S(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.rk(1)*sim.K(1) - sim.PiQ(1) ...
- (sim.Y(1) - sim.w(1)*sim.L(1) - sim.rk(1)*sim.K(1)) - sim.w(1)*sim.L(1) + sim.Chi(1)

sim.S(1) + sim.Xi(1) + sim.C(1) + sim.Q(1) * sim.IH(1) - sim.PiQ(1) - sim.Y(1) + sim.Chi(1)

sim.Y(1) - (sim.C(1) + sim.Q(1) * sim.IH(1) - sim.PiQ(1) + sim.Chi(1) + sim.Xi(1))

=> 
sim.Q(1) * sim.IH(1) - sim.PiQ(1) - sim.I(1) % not correct!

sim.Q(1) * sim.IH(1) - sim.PiQ(1) - sim.I(1) 

sim.PiQ(1) = (sim.Q(1) - 1).*param.solve_for_iota(sim.Q(1)).*sim.K(1) - param.Phi(param.solve_for_iota(sim.Q(1))).*sim.K(1)

sim.I(1) = param.solve_for_iota(sim.Q(1)).*sim.K(1) + param.Phi(param.solve_for_iota(sim.Q(1))).*sim.K(1)

SO: 
sim.PiQ(1) = sim.Q(1)*param.solve_for_iota(sim.Q(1))*sim.K(1) - sim.I(1)

Therefore: We need 
IH = param.solve_for_iota(Q)*K 
and that's not the case

%}


