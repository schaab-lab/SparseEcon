function [sim, G, G_dense] = transition(x, t, z, ss, G, G_dense, param, query)


%% AGGREGATE TRANSITION PATH
sim = macro_block(x, t, z, ss, param);

if nargin > 7 && ~any(ismember(query, {'u', 'c', 's', 'm', 'V', 'g', 'markets', 'all'}))
    sim = sim.(query);
    return; 
end


%% PREALLOCATE
[sim.V, sim.A, sim.g, sim.c, sim.s, sim.m] = deal(cell(param.N, 1));
[sim.C, sim.B, sim.S, sim.KH, sim.MH, sim.IH, sim.Chi, sim.Xi] = deal(zeros(param.N, 1));


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
       
    % POLICY FUNCTIONS
    ltau  = 15;
    ltau0 = sim.rk(n) * (param.kmax*0.999)^(1-ltau);
    sim.xi{n} = param.xi * ltau0 * G.k .^ ltau;

    G.income_a = sim.r(n) * G.a + sim.rk(n) * G.k - sim.xi{n} ...
                 + (1-param.tau_lab) * sim.w(n) .* param.zz .* sim.N(n) + sim.tau(n);
    G.income_k = - param.delta * G.k;
    G.N = sim.N(n); G.piw = sim.piw(n);
    
    hjb = HJB(V, G, param);
    
    % FD OPERATORS
    Asc = cell(param.discrete_types, 1); Asi = Asc; Ami = Asc; Amk = Asc;
    for j = 1:param.discrete_types
        Asc{j} = FD_operator(G, hjb.sc(:, j),   zeros(G.J, 1), 1, num2str(j));
        Asi{j} = FD_operator(G, hjb.si(:, j),   zeros(G.J, 1), 1, num2str(j));
        Ami{j} = FD_operator(G, hjb.iota(:, j), zeros(G.J, 1), 2, num2str(j));
        Amk{j} = FD_operator(G, G.income_k,     zeros(G.J, 1), 2, num2str(j));
    end

    A = blkdiag(Asc{1} + Asi{1} + Ami{1} + Amk{1}, Asc{2} + Asi{2} + Ami{2} + Amk{2}) + Az;

    B = (1/param.dt + sim.rho(n))*speye(2*G.J) - A;
    b = hjb.u(:) + V(:) / param.dt;

    % SOLVE LINEAR SYSTEM
    V_new = B \ b;
    V = [V_new(1:G.J), V_new(1+G.J:end)];
    if ~isreal(V), disp('Complex values detected!'); sim = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n} = V; sim.u{n} = hjb.u; sim.c{n} = hjb.c; sim.m{n} = hjb.m; sim.s{n} = hjb.s; 
    sim.iota{n} = hjb.iota; sim.sc{n} = hjb.sc; sim.si{n} = hjb.si; sim.A{n} = A;
    
end


%% SOLVE KF FORWARDS
% dgdtA = [Aa1' * g_t{n}(:, 1); Aa2' * g_t{n}(:, 2)] + ...
%         [Ak1' * g_t{n}(:, 1); Ak2' * g_t{n}(:, 2)] + ...
%          Az' * [g_t{n}(:, 1); g_t{n}(:, 2)];
% dgdtA = [dgdtA(1:GDense.J), dgdtA(GDense.J+1:end)];
Az_dense = [-speye(G_dense.J)*param.la1,  speye(G_dense.J)*param.la1; ...
             speye(G_dense.J)*param.la2, -speye(G_dense.J)*param.la2];

% Birth process:
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
        B = 1/param.dt * speye(param.discrete_types*G_dense.J) - AT;
        b = sim.g{n}(:) / param.dt;

        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    else
        % Explicit:
        gg = sim.g{n}(:) + param.dt * AT * sim.g{n}(:);
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    end
    
    % Ensure positive density:
    sim.g{n+1}(sim.g{n+1}(:, 1) < 0, 1) = 0; 
    sim.g{n+1}(sim.g{n+1}(:, 2) < 0, 2) = 0; 
    sim.g{n+1}(:, 1) = ss.mass(1) * sim.g{n+1}(:, 1) / sum(sum(sim.g{n+1}(:, 1) .* G_dense.dx));
    sim.g{n+1}(:, 2) = ss.mass(2) * sim.g{n+1}(:, 2) / sum(sum(sim.g{n+1}(:, 2) .* G_dense.dx));
    sim.mass(n+1, :) = sum( sim.g{n+1} .* G_dense.dx );
    
    assert(all(all( sim.g{n+1} >= 0 )));
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - 1) > 1e-4
        error('Sim KF not mass-preserving.\n'); 
    end
    
end

if nargin > 7 && ~any(ismember(query, {'markets', 'all'}))
    sim = sim.(query);
    return; 
end


%% AGGREGATION
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
    u1z_dense  = param.zz .* param.u1(c_dense);
    
    sim.B(n)  = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.KH(n) = sum(sum(G_dense.k .* sim.g{n} .* G_dense.dx));
    sim.C(n)  = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n)  = sum(sum(s_dense .* sim.g{n} .* G_dense.dx));
    sim.S2(n) = sum(sum((sc_dense+si_dense) .* sim.g{n} .* G_dense.dx));
    sim.M1(n) = sum(sum(m_dense .* sim.g{n} .* G_dense.dx));
    sim.M2(n) = sum(sum((mi_dense+mk_dense) .* sim.g{n} .* G_dense.dx));
    sim.IH(n) = sum(sum(iota_dense .* sim.g{n} .* G_dense.dx));
    sim.Chi(n)= sum(sum(adjcostfn(iota_dense, G_dense.k, param) .* sim.g{n} .* G_dense.dx));
    sim.Xi(n) = sum(sum(xi_dense .* sim.g{n} .* G_dense.dx)); 
    sim.MH(n) = sum(sum(u1z_dense .* sim.g{n} .* G_dense.dx));
    
end

% sim.K2(1) = ss.K; sim.K3(1) = ss.K; sim.B2(1) = ss.B; sim.B3(1) = ss.B;
% for n = 1:param.N-1
%     sim.K2(n+1) = sim.K2(n) + param.dt * sim.M(n);
%     sim.K3(n+1) = sim.K3(n) + param.dt * (sim.M3(n) - param.delta *sim.K3(n));
%     sim.B2(n+1) = sim.B2(n) + param.dt * sim.S(n);
%     sim.B3(n+1) = sim.B3(n) + param.dt * sim.S2(n);
% end


%% MARKET CLEARING
sim.diff_K = sim.KH - sim.K;
sim.diff_I = sim.IH - sim.I;
sim.diff_S = sim.S;
sim.diff_B = sim.B - param.gov_bond_supply;
sim.diff_Y = sim.Y - sim.C - sim.I - sim.Chi - sim.Xi - sim.G;
sim.diff_M = sim.M - sim.MH;

% sim.excess_MPL = sim.N - ( (1-param.alpha) * sim.Y ./ sim.w );
% sim.excess_MPK = sim.K - ( param.alpha * sim.Y ./ sim.rk );

sim.diff_markets = [sim.diff_Y, sim.diff_B, sim.diff_M, sim.diff_S, sim.diff_K, sim.diff_I];

if nargin > 7 && any(ismember(query, {'markets'}))
    % sim = sim.(query);
    % Cannot use diff_B or diff_K because at n=0, they are =0 by ss init!
    sim = [sim.diff_S; sim.diff_I; sim.diff_M];
end

end




