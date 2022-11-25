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



