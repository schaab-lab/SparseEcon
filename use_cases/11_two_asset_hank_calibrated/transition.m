function [sim, G, G_dense] = transition(x, z, ss, G, G_dense, param, query)


%% AGGREGATE TRANSITION PATH
sim = macro_block(x, z, ss, param);

if nargin > 6 && ~any(ismember(query, {'u', 'c', 's', 'm', 'V', 'g', 'markets', 'all'}))
    sim = sim.(query); return; 
end


%% PREALLOCATE
[sim.V, sim.A, sim.g, sim.c, sim.s, sim.m, sim.iota, sim.sc, sim.si, sim.xi] = deal(cell(param.N, 1));
[sim.C, sim.B, sim.S, sim.KH, sim.IH, sim.Lambda, sim.Chi, sim.Xi, ...
    sim.LH, sim.M, sim.I2, sim.B2, sim.B3, sim.S2, sim.M2, sim.K2, sim.K3] = deal(zeros(param.N, 1));


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
       
    % Policy functions:
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
    
    % Finite-difference operators:
    Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
    for j = 1:param.discrete_types
        Asc{j} = FD_operator(G, hjb.sc(:, j),   zeros(G.J, 1), 1, num2str(j));
        Asi{j} = FD_operator(G, hjb.si(:, j),   zeros(G.J, 1), 1, num2str(j));
        Ami{j} = FD_operator(G, hjb.iota(:, j), zeros(G.J, 1), 2, num2str(j));
        Amk{j} = FD_operator(G, G.income_k,    zeros(G.J, 1), 2, num2str(j));
    end

    A = blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az;

    B = (1/param.dt(n) + param.rho + param.deathrate)*speye(2*G.J) - A;
    b = hjb.u(:) + V(:) / param.dt(n);

    % Solve linear system:
    V_new = B \ b;
    V = [V_new(1:G.J), V_new(1+G.J:end)];    
    if ~isreal(V), disp('Complex values detected!'); DIFF = NaN(1); return; end
    
    % Record data:
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
        B = (1/param.dt(n) + param.deathrate) * speye(param.discrete_types*G_dense.J) - AT;
        b = sim.g{n}(:) / param.dt(n) + param.deathrate * birth_ID(:) / G_dense.dx;

        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    else
        if n < param.N
            t_KF  = linspace(param.t(n), param.t(n+1), param.reso_sim_KF+1);
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
    
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - sum(sum(sim.g{n} * G_dense.dx ))) > 1e-6
        fprintf('KF not preserving mass.\n');
    end    
    
end

if nargin > 6 && ~any(ismember(query, {'markets', 'all'}))
    sim = sim.(query);
    return; 
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
    sim.B(n)       = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.KH(n)      = sum(sum(G_dense.k .* sim.g{n} .* G_dense.dx));
    sim.C(n)       = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n)       = sum(sum(s_dense .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.B(n);
    sim.S2(n)      = sum(sum((sc_dense+si_dense) .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.B(n);
    sim.M(n)       = sum(sum(m_dense .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.KH(n);
    sim.M2(n)      = sum(sum((mi_dense+mk_dense) .* sim.g{n} .* G_dense.dx)) - param.deathrate * sim.KH(n);
    sim.IH(n)      = sum(sum(iota_dense .* sim.g{n} .* G_dense.dx));
    sim.Chi(n)     = sum(sum(adjcostfn(iota_dense, G_dense.k, param) .* sim.g{n} .* G_dense.dx));
    sim.Xi(n)      = sum(sum(xi_dense .* sim.g{n} .* G_dense.dx)); 
    sim.Lambda(n)  = sum(sum( (param.zz .* param.u1(c_dense)) .* sim.g{n} .* G_dense.dx));
    
end
sim.iotaQ = param.solve_for_iota(sim.Q);
sim.Phi = param.Phi(sim.iotaQ) .* sim.K;
sim.IQ = sim.iotaQ .* sim.K;
sim.I2 = param.gross_total_investment_expenditure(sim.Q, sim.K, 0);
sim.Q2 = 1 + param.Phi_prime(sim.dK ./ sim.K + param.delta);
sim.Q3 = 1 + param.Phi_prime(sim.iotaQ);
sim.M3 = sim.gross_total_capital_accumulation;

sim.K2(1) = ss.K; sim.K3(1) = ss.K; sim.B2(1) = ss.B; sim.B3(1) = ss.B;
for n = 1:param.N-1
    sim.K2(n+1) = sim.K2(n) + param.dt(n) * sim.M(n);
    sim.K3(n+1) = sim.K3(n) + param.dt(n) * (sim.M3(n) - param.delta *sim.K3(n));
    sim.B2(n+1) = sim.B2(n) + param.dt(n) * sim.S(n);
    sim.B3(n+1) = sim.B3(n) + param.dt(n) * sim.S2(n);
end


%% GAPS
sim.diff_K = sim.KH - sim.K;
sim.diff_I = sim.gross_total_capital_accumulation - sim.IH;
sim.diff_S = sim.S;
sim.diff_B = sim.B - param.gov_bond_supply;
sim.diff_Y = sim.Y - sim.C - sim.I - sim.Chi - sim.Xi - sim.G;
sim.diff_L = param.v1(sim.H) - (param.epsilonW-1)/param.epsilonW*(1-param.tau_lab)*sim.w.*sim.Lambda; %sim.M - sim.MH;
sim.diff_T = sim.tau - (sim.Pi + param.tau_lab * sim.w .* sim.L ...
             - param.UI .* sim.U - sim.G - sim.r*param.gov_bond_supply);

sim.diff_MPL = sim.L - ( (1-param.alpha) * sim.mc .* sim.Y ./ ((1-param.tau_empl)*sim.w) );
sim.diff_MPK = sim.K - ( param.alpha * sim.mc .* sim.Y ./ sim.rk );

sim.diff_markets = [sim.diff_Y, sim.diff_B, sim.diff_L, sim.diff_S, sim.diff_K, sim.diff_I, ...
                    sim.diff_T, sim.diff_MPL, sim.diff_MPK];

if nargin > 6 && any(ismember(query, {'markets'}))
    % Careful with diff_B or diff_K: at n=0, they are =0 by ss initialization
    sim = [interp1(param.t, sim.diff_S, param.nodes)'; ...
           interp1(param.t, sim.diff_K, param.nodes)'; ... %sim.diff_I; ...
           interp1(param.t, sim.diff_L, param.nodes)'];
end


end




