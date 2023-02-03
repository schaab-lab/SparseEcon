function [sim, G, G_dense] = transition(x, z, ss, G, G_dense, param, query)


%% MACRO BLOCK: PRE
sim = macro_block_pre(x, z, ss, param);

if nargin > 6 && ~any(ismember(query, {'u', 'c', 's', 'm', 'V', 'g', 'adjoint', 'fake_news', 'markets', 'all'}))
    sim = sim.(query);
    return;
end


%% PREALLOCATE
[sim.V, sim.A, sim.g, sim.c, sim.s, sim.m] = deal(cell(param.N, 1));
[sim.C, sim.B, sim.S, sim.M] = deal(zeros(param.N, 1));


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
    
    G.income = sim.r(n) * G.a + sim.w(n) .* param.zz * sim.N(n);
    G.N = sim.N(n);
    
    % POLICY FUNCTIONS
    num0 = 1e-8; % numerical 0 for upwind scheme
    
    [VaF, VaB] = deal(zeros(G.J, param.discrete_types));
    for j = 1:param.discrete_types
        VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    
    VaF(G.grid(:, 1) == 1, :) = param.u1(G.income(G.grid(:, 1) == 1, :));
    VaB(G.grid(:, 1) == 0, :) = param.u1(G.income(G.grid(:, 1) == 0, :));
    
    cF = param.u1inv(VaF);
    cB = param.u1inv(VaB);
    c0 = G.income;
    
    sF = G.income - cF;
    sB = G.income - cB;
    
    IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
    IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
    I0 = ~IF & ~IB;
    
    s = sF.*IF + sB.*IB;
    c = cF.*IF + cB.*IB + c0.*I0;
    u = param.u(c) - param.v(G.N);
    m = param.zz .* param.u1(c);
    
    % CONSTRUCT FD OPERATORS
    [Aa1, const1] = FD_operator(G, s(:, 1), zeros(G.J, 1), 1, '1');
    [Aa2, const2] = FD_operator(G, s(:, 2), zeros(G.J, 1), 1, '2');
    
    A = blkdiag(Aa1, Aa2) + Az;
    const = [const1; const2];
    
    B = (1/param.dt + sim.rho(n))*speye(param.discrete_types*G.J) - A;
    b = u(:) + V(:) / param.dt + const;    
    
    % SOLVE LINEAR SYSTEM
    V_new = B\b;
    V = [V_new(1:G.J), V_new(1+G.J:end)];
    if ~isreal(V), disp('Complex values detected!'); sim = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n} = V; sim.u{n} = u; sim.c{n} = c; sim.s{n} = s; sim.A{n} = A; sim.m{n} = m;
    
end

if nargin > 6 && ~any(ismember(query, {'g', 'adjoint', 'markets', 'all'}))
    if ~any(ismember(query, {'fake_news'})), sim = sim.(query); end
    return;
end


%% SOLVE KF FORWARDS
Az_dense = [-speye(G_dense.J)*param.la1,  speye(G_dense.J)*param.la1; ...
             speye(G_dense.J)*param.la2, -speye(G_dense.J)*param.la2];

for n = 1:param.N

    Aa_dense1 = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 1), zeros(G_dense.J, 1), 1, '1');
    Aa_dense2 = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 2), zeros(G_dense.J, 1), 1, '2');
    
    AT = (blkdiag(Aa_dense1, Aa_dense2) + Az_dense)';
    if isequal(query, 'adjoint')
        sim = speye(G_dense.J*param.discrete_types) + param.dt * AT;
        return;
    end
    
    if param.implicit_g
        % Implicit:
        B = 1/param.dt * speye(G_dense.J*param.discrete_types) - AT;
        b = sim.g{n}(:) / param.dt;
        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    else
        % Explicit:
        gg = sim.g{n}(:) + param.dt * AT * sim.g{n}(:);
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    end
    
    % Ensure positive density:
    sim.g{n+1}(sim.g{n+1} < 0) = 0;    
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - sum(sum(sim.g{n} * G_dense.dx ))) > 1e-6
        fprintf('KF not preserving mass.\n');
    end    
    
end

if nargin > 6 && ~any(ismember(query, {'markets', 'all'}))
    sim = sim.(query);
    return; 
end


%% AGGREGATION
for n = 1:param.N
    c_dense = G.BH_dense * sim.c{n};
    s_dense = G.BH_dense * sim.s{n};
    m_dense = G.BH_dense * sim.m{n};
    
    sim.mass(n, :) = sum(sim.g{n}*G_dense.dx);
    
    sim.B(n) = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.C(n) = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n) = sum(sum(s_dense .* sim.g{n} .* G_dense.dx));
    
    sim.M(n) = sum(sum(m_dense .* sim.g{n} .* G_dense.dx));
    % for j = 1:numel(aggregates)
    %     sim.([aggregates{j}])(n) = ...
    %         sum(sum(eval([lower(aggregates{j}), '_dense']) .* sim.g{n} .* G_dense.dx));
    % end
end


%% MACRO BLOCK: POST
c = [sim.B; sim.C; sim.S; sim.M];
aggregates = {'B', 'C', 'S', 'M'};
sim = macro_block_post(x, c, z, sim, ss, param, aggregates);


%% OUTPUT
sim.excess_bonds = sim.B;
sim.excess_goods = sim.Y - sim.C;
sim.excess_union = sim.MX - sim.M;
sim.excess_saving = sim.S;

sim.diff_markets = [sim.excess_goods, sim.excess_bonds, sim.excess_union, sim.excess_saving];

if nargin > 6 && ~any(ismember(query, {'all'}))
    sim = [sim.excess_saving; sim.excess_union];
end


end




