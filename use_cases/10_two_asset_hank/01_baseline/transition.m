function [DIFF, G, G_dense, sim] = transition(PHI, G, G_dense, shock, ss, param)

%% SETUP
sim.N  = param.N;
sim.t  = param.t;
sim.dt = param.dt;

% PREALLOCATE
sim.V  = cell(param.N, 1);
sim.A  = cell(param.N, 1);
sim.g  = cell(param.N, 1);
sim.c  = cell(param.N, 1);
sim.s  = cell(param.N, 1);
sim.m  = cell(param.N, 1);
sim.sc = cell(param.N, 1);
sim.si = cell(param.N, 1);
sim.xi = cell(param.N, 1);
sim.iota = cell(param.N, 1);

sim.KH = zeros(param.N, 1);
sim.LH = zeros(param.N, 1);
sim.C  = zeros(param.N, 1);
sim.S  = zeros(param.N, 1);
sim.M  = zeros(param.N, 1);
sim.B  = zeros(param.N, 1);
sim.I  = zeros(param.N, 1);
sim.IH = zeros(param.N, 1);
sim.Xi = zeros(param.N, 1);
sim.B2 = zeros(param.N, 1);
sim.B3 = zeros(param.N, 1);
sim.S2 = zeros(param.N, 1);
sim.M2 = zeros(param.N, 1);
sim.K2 = zeros(param.N, 1);
sim.K3 = zeros(param.N, 1);
sim.Chi = zeros(param.N, 1);
sim.Lambda = zeros(param.N, 1);
sim.Phi = zeros(param.N, 1);
sim.iotaQ = zeros(param.N, 1);
sim.IQ = zeros(param.N, 1);

sim.Q2 = zeros(param.N, 1);
sim.Q3 = zeros(param.N, 1);


%% AGGREGATE TRANSITION PATH
X = basis_fun_irf([], reshape(PHI, [1, numel(PHI)]), param.H(1), param.H(2), param.bfun_type, sim.t, "get_function");

sim.i  = X(:, 1);
sim.K  = X(:, 2);
sim.L  = X(:, 3);

sim.Z   = zeros(param.N, 1);
sim.eps = zeros(param.N, 1);
sim.G   = zeros(param.N, 1);
switch param.shock_type
    case 'monetary'
        sim.eps = shock;
        
    case 'productivity'
        sim.Z  = shock;
        
    case 'fiscal'
        sim.G  = shock;
end

sim.Y  = exp(sim.Z) .* sim.K.^param.alpha .* sim.L.^(1-param.alpha);

sim.pi = (sim.i - sim.eps - ss.r) / param.lambda_pi;
sim.r  = sim.i - sim.pi;

sim.P   = ones(param.N, 1);
sim.dY  = zeros(param.N, 1);
sim.dpi = zeros(param.N, 1);
for n = 1:param.N-1
    sim.P(n+1) = sim.P(n) .* (1+sim.dt(n)*sim.pi(n));
    sim.dY(n)  = (sim.Y(n+1) - sim.Y(n)) / sim.dt(n);
    sim.dpi(n) = (sim.pi(n+1) - sim.pi(n)) / sim.dt(n);
end

%%% Using the nominal riskfree rate as the firm's discount rate implies
%%% that we don't have to guess 'ik' beforehand. 
sim.mc = (param.epsilonF-1)/param.epsilonF * (1 + param.chiF/(param.epsilonF-1) ...
         .* (sim.pi.*(sim.i - sim.pi - sim.dY./sim.Y) - sim.dpi));
sim.rk = param.alpha * sim.mc .* sim.Y ./ sim.K;
sim.ik = sim.rk .* sim.P;
sim.w  = 1/(1-param.tau_empl) * ( exp(sim.Z) .* sim.mc * (param.alpha^param.alpha * (1-param.alpha)^(1-param.alpha)) ...
         .* (sim.rk).^(-param.alpha) ).^(1./(1-param.alpha));

sim.Pi = (1-sim.mc) .* sim.Y;
% sim.L  = sim.Y .* (param.alpha/(1-param.alpha) .* (1-param.tau_empl) .* sim.w ./ sim.rk).^(-param.alpha);
% sim.K  = param.alpha./(1-param.alpha) .* (1-param.tau_empl) .* sim.L .* sim.w ./ sim.rk;
% sim.L = (1-param.alpha) * sim.mc .* sim.Y ./ ((1-param.tau_empl)*sim.w);
% sim.K = param.alpha * sim.mc .* sim.Y ./ sim.rk;

sim.dK = zeros(param.N, 1);
for n = 1:param.N-1
    sim.dK(n)  = (sim.K(n+1) - sim.K(n)) / sim.dt(n);
end

sim.gross_total_capital_accumulation = sim.dK + param.delta * sim.K;
sim.Q = param.solve_for_Q_from_cap_accumulation(sim.gross_total_capital_accumulation, sim.K);

sim.gross_total_investment_expenditure = param.gross_total_investment_expenditure(sim.Q, sim.K, 0);
sim.iotaQ = param.solve_for_iota(sim.Q);
% assert( max(abs( sim.gross_total_capital_accumulation - param.gross_total_capital_accumulation(sim.Q, sim.K, 0) )) < 1e-8 );
% assert( max(abs( sim.Q - param.solveForQ(sim.iotaQ) )) < 1e-8 );
% assert( max(abs( sim.iotaQ - sim.gross_total_capital_accumulation ./ sim.K )) < 1e-8 );

sim.I   = sim.gross_total_investment_expenditure;
sim.PiQ = param.PiQ(sim.Q, sim.K, zeros(param.N, 1));

sim.U = param.la2 / (param.la1 + param.la2);
sim.H = sim.L ./ (1-sim.U);

sim.tau = sim.Pi + param.tau_lab * sim.w .* sim.L - param.UI * sim.U - sim.G - sim.r * param.gov_bond_supply;


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
       
    % POLICY FUNCTIONS
    cap_income_liquid   = sim.rk(n) + sim.PiQ(n) / sim.K(n);
    cap_income_illiquid = param.deathrate - param.delta;

    ltau  = 15;
    ltau0 = cap_income_liquid * (param.kmax*0.999)^(1-ltau);
    sim.xi{n} = param.xi * ltau0 * G.k .^ ltau;

    G.income_a = (sim.r(n) + param.deathrate) .* G.a ...
                 + cap_income_liquid .* G.k - sim.xi{n} ...
                 + (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.H(n) + sim.tau(n) + param.UI .* [1, 0];
    G.income_k = cap_income_illiquid .* G.k;
    G.Q = sim.Q(n); G.H = sim.H(n);

    hjb = HJB(V, G, param);
    
    c_diff(n) = max(max(abs(hjb.c - ss.c)));
    
    % FD OPERATORS
    Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
    for j = 1:param.discrete_types
        Asc{j} = FD_operator(G, hjb.sc(:, j),   zeros(G.J, 1), 1, num2str(j));
        Asi{j} = FD_operator(G, hjb.si(:, j),   zeros(G.J, 1), 1, num2str(j));
        Ami{j} = FD_operator(G, hjb.iota(:, j), zeros(G.J, 1), 2, num2str(j));
        Amk{j} = FD_operator(G, G.income_k,    zeros(G.J, 1), 2, num2str(j));
    end

    A = blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az;

    B = (1/sim.dt(n) + param.rho + param.deathrate)*speye(2*G.J) - A;
    b = [hjb.u(:, 1); hjb.u(:, 2)] + [V(:, 1); V(:, 2)] / sim.dt(n);

    % SOLVE LINEAR SYSTEM
    V_new = B\b; V_diff(n) = max(max(abs([V_new(1:G.J), V_new(1+G.J:end)]-V)));

    V = [V_new(1:G.J), V_new(1+G.J:end)];
    
    if ~isreal(V), disp('Complex values detected!'); DIFF = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n} = V; sim.u{n} = hjb.u; sim.c{n} = hjb.c; sim.m{n} = hjb.m; sim.s{n} = hjb.s; 
    sim.iota{n} = hjb.iota; sim.sc{n} = hjb.sc; sim.si{n} = hjb.si;
end


%% SOLVE KF FORWARDS
% dgdtA = [Aa1' * g_t{n}(:, 1); Aa2' * g_t{n}(:, 2)] + ...
%         [Ak1' * g_t{n}(:, 1); Ak2' * g_t{n}(:, 2)] + ...
%          Az' * [g_t{n}(:, 1); g_t{n}(:, 2)];
% dgdtA = [dgdtA(1:GDense.J), dgdtA(GDense.J+1:end)];
Az_dense = [-speye(G_dense.J)*param.la1,  speye(G_dense.J)*param.la1; ...
             speye(G_dense.J)*param.la2, -speye(G_dense.J)*param.la2];

% Birth process:
[~, birth_idx] = min( (G_dense.a).^2 + (G_dense.k).^2 );
if G_dense.a(birth_idx, :) < 0, birth_idx2 = birth_idx+1; elseif G_dense.a(birth_idx, :) > 0, birth_idx2 = birth_idx-1; end
birth_ID = zeros(G_dense.J, 2); 
if G_dense.a(birth_idx) == 0
    birth_ID(birth_idx, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)];
else
    birth_ID(birth_idx, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)] ...
                            * ( abs(G_dense.a(birth_idx2)) / ...
                              ( abs(G_dense.a(birth_idx)) + abs(G_dense.a(birth_idx2))));
    birth_ID(birth_idx2, :) = [param.la2/(param.la1+param.la2), param.la1/(param.la1+param.la2)] ...
                            * ( abs(G_dense.a(birth_idx)) / ...
                              ( abs(G_dense.a(birth_idx)) + abs(G_dense.a(birth_idx2))));
end

for n = 1:param.N
    
    Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
    for j = 1:param.discrete_types
        Asc{j} = FD_operator(G_dense, G.BH_dense * sim.sc{n}(:, j),   zeros(G_dense.J, 1), 1, num2str(j));
        Asi{j} = FD_operator(G_dense, G.BH_dense * sim.si{n}(:, j),   zeros(G_dense.J, 1), 1, num2str(j));
        Ami{j} = FD_operator(G_dense, G.BH_dense * sim.iota{n}(:, j), zeros(G_dense.J, 1), 2, num2str(j));
        Amk{j} = FD_operator(G_dense, G.BH_dense * G.income_k,        zeros(G_dense.J, 1), 2, num2str(j));
    end
    AT = (blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az_dense)';
    
    if param.implicit_g
        % Implicit:
        B = (1/sim.dt(n) + param.deathrate) * speye(param.discrete_types*G_dense.J) - AT;
        b = [sim.g{n}(:, 1); sim.g{n}(:, 2)]/sim.dt(n) + param.deathrate * [birth_ID(:, 1); birth_ID(:, 2)] / G_dense.dx;

        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    else
        if n < param.N
            t_KF  = linspace(sim.t(n), sim.t(n+1), param.reso_sim_KF+1);
            dt_KF = t_KF(2) - t_KF(1);
        end
        g_KF    = cell(1, param.reso_sim_KF+1);
        g_KF{1} = sim.g{n};

        for m = 1:param.reso_sim_KF
            dgdtA = [Asc{1}' * g_KF{m}(:, 1); Asc{2}' * g_KF{m}(:, 2)] + ...
                    [Asi{1}' * g_KF{m}(:, 1); Asi{2}' * g_KF{m}(:, 2)] + ...
                    [Ami{1}' * g_KF{m}(:, 1); Ami{2}' * g_KF{m}(:, 2)] + ...
                    [Amk{1}' * g_KF{m}(:, 1); Amk{2}' * g_KF{m}(:, 2)] + ...
                     Az_dense' * [g_KF{m}(:, 1); g_KF{m}(:, 2)];
            % dgdtA = [Aa1' * gKF{m}(:, 1); Aa2' * gKF{m}(:, 2)] + ...
            %         [Ak1' * gKF{m}(:, 1); Ak2' * gKF{m}(:, 2)] + ...
            %         [At'  * gKF{m}(:, 1); At'  * gKF{m}(:, 2)] + ...
            %          Az' * [gKF{m}(:, 1); gKF{m}(:, 2)];
            dgdtA = [dgdtA(1:G_dense.J), dgdtA(G_dense.J+1:end)];

            g_KF{m+1} = g_KF{m} + dgdtA * dt_KF ...
                       + param.deathrate * birth_ID/G_dense.dx * dt_KF ...
                       - param.deathrate * g_KF{m} * dt_KF;
        end
        sim.g{n+1} = g_KF{param.reso_sim_KF+1};
    end
    
    ID1 = sim.g{n+1}(:, 1)<0;
    ID2 = sim.g{n+1}(:, 2)<0;
    sim.g{n+1}(ID1, 1) = 0; 
    sim.g{n+1}(ID2, 2) = 0; 
    sim.g{n+1}(:, 1) = ss.mass(1) * sim.g{n+1}(:, 1) / sum(sum(sim.g{n+1}(:, 1) .* G_dense.dx));
    sim.g{n+1}(:, 2) = ss.mass(2) * sim.g{n+1}(:, 2) / sum(sum(sim.g{n+1}(:, 2) .* G_dense.dx));
    sim.mass(n+1, :) = sum( sim.g{n+1} .* G_dense.dx );
    
    assert(all(all( sim.g{n+1} >= 0 )));   
    if abs( sum(sum(sim.g{n+1} .* G_dense.dx )) - 1) > 1e-7, error('Sim KF not mass-preserving.\n'); end
    
    % g_diff(n) = max(max(abs( sim.g{n} - ss.g )));    
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - sum(sum(sim.g{n} * G_dense.dx ))) > 1e-6
        fprintf('KF not preserving mass.\n');
    end    
    
end


%% AGGREGATION & MARKET CLEARING
for n = 1:param.N    
    
    c_dense    = G.BH_dense * sim.c{n};
    s_dense    = G.BH_dense * sim.s{n};
    si_dense   = G.BH_dense * sim.si{n};
    sc_dense   = G.BH_dense * sim.sc{n};
    m_dense    = G.BH_dense * sim.m{n};
    xi_dense   = G.BH_dense * sim.xi{n};
    iota_dense = G.BH_dense * sim.iota{n};
    mi_dense   = iota_dense;
    mk_dense   = G.BH_dense * G.income_k;
    
    sim.mass(n, :) = sum(sim.g{n}*G_dense.dx);    
    sim.B(n)      = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.KH(n)     = sum(sum(G_dense.k .* sim.g{n} .* G_dense.dx));
    sim.C(n)      = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n)      = sum(sum(s_dense .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.B(n);
    sim.S2(n)     = sum(sum((sc_dense+si_dense) .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.B(n);
    sim.M(n)      = sum(sum(m_dense .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.KH(n);
    sim.M2(n)     = sum(sum((mi_dense+mk_dense) .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.KH(n);
    sim.I(n)      = param.gross_total_investment_expenditure(sim.Q(n), sim.K(n), 0);
    sim.IH(n)     = sum(sum(iota_dense .* sim.g{n} .* G_dense.dx));
    sim.Chi(n)    = sum(sum(adjcostfn(iota_dense, G_dense.k, param) .* sim.g{n} .* G_dense.dx));
    sim.Xi(n)     = sum(sum(xi_dense .* sim.g{n} .* G_dense.dx)); 
    sim.Lambda(n) = sum(sum( (param.zz .* param.u1(c_dense)) .* sim.g{n} .* G_dense.dx));
    
end
sim.iotaQ = param.solve_for_iota(sim.Q);
sim.Phi = param.Phi(sim.iotaQ) .* sim.K;
sim.IQ = sim.iotaQ .* sim.K;

sim.Q2 = 1 + param.Phi_prime(sim.dK ./ sim.K + param.delta);
sim.Q3 = 1 + param.Phi_prime(sim.iotaQ);


sim.excess_capital = sim.KH - sim.K;
sim.excess_saving  = sim.S;
sim.excess_bonds   = sim.B - param.gov_bond_supply;
sim.excess_labor   = param.v1(sim.H) - (param.epsilonW-1)/param.epsilonW*(1-param.tau_lab)*sim.w.*sim.Lambda;
sim.excess_goods   = sim.Y - sim.C - sim.I - sim.Chi - sim.Xi - sim.G;
% sim.excess_tax     = sim.tau - (sim.Pi + param.tau_lab * sim.w .* sim.L ...
%                                 - param.UI .* sim.U - sim.G - sim.r*param.gov_bond_supply);
sim.excess_cap_production = sim.gross_total_capital_accumulation - sim.IH;

sim.excess_MPL = sim.L - ( (1-param.alpha) * sim.mc .* sim.Y ./ ((1-param.tau_empl)*sim.w) );
sim.excess_MPK = sim.K - ( param.alpha * sim.mc .* sim.Y ./ sim.rk );

sim.M3 = sim.gross_total_capital_accumulation;

sim.K2(1) = ss.K; sim.K3(1) = ss.K; sim.B2(1) = ss.B; sim.B3(1) = ss.B;
for n = 1:param.N-1
    sim.K2(n+1) = sim.K2(n) + sim.dt(n) * sim.M(n);
    sim.K3(n+1) = sim.K3(n) + sim.dt(n) * (sim.M3(n) - param.delta *sim.K3(n));
    sim.B2(n+1) = sim.B2(n) + sim.dt(n) * sim.S(n);
    sim.B3(n+1) = sim.B3(n) + sim.dt(n) * sim.S2(n);
end



%% COLLOCATION POINTS

DIFF_i  = interp1(sim.t, sim.excess_saving,   param.nodes);
DIFF_K  = interp1(sim.t, sim.excess_capital, param.nodes);
DIFF_L  = interp1(sim.t, sim.excess_labor,   param.nodes);

DIFF = [DIFF_i, DIFF_K, DIFF_L]';


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

end




