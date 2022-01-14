function [diff, G, G_dense, sim] = transition(PHI, G, G_dense, shock, ss, param)

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
sim.C = zeros(param.N, 1);
sim.B = zeros(param.N, 1);
sim.S = zeros(param.N, 1);


%% AGGREGATE TRANSITION PATH
X = basis_fun_irf([], reshape(PHI, [1, numel(PHI)]), param.H(1), param.H(2), ...
    param.bfun_type, sim.t, "get_function");

sim.r = X(:,1); 
sim.Z = shock;

sim.L = param.L * ones(param.N, 1);
sim.Y = exp(sim.Z) .* sim.L;
sim.w = exp(sim.Z);


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

[Az, const_z] = FD_operator(G, param.theta_z * (param.zmean - G.z), param.sig_z*ones(G.J,1), 2);

for n = param.N:-1:1
    
    G.income = sim.r(n) * G.a + sim.w(n) .* G.z;
    
    % POLICY FUNCTIONS
    num0 = 1e-8; % numerical 0 for upwind scheme

    VaF = deriv_sparse(G, V, 1, 'D1F');
    VaB = deriv_sparse(G, V, 1, 'D1B');
    
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
    u = param.u(c);
    
    % CONSTRUCT FD OPERATORS
    [Aa, const_a] = FD_operator(G, s, zeros(G.J,1), 1);
    
    A = Aa + Az;
    const = const_a + const_z;

    B = (1/sim.dt(n) + param.rho)*speye(G.J) - A;
    b = u + V / sim.dt(n) + const;    
    
    % SOLVE LINEAR SYSTEM
    V = B\b;    
    if ~isreal(V), disp('Complex values detected!'); diff = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n} = V; sim.u{n} = u; sim.c{n} = c; sim.s{n} = s;
    
end


%% SOLVE KF FORWARDS
% dgdtA = [Aa1' * g_t{n}(:, 1); Aa2' * g_t{n}(:, 2)] + ...
%         [Ak1' * g_t{n}(:, 1); Ak2' * g_t{n}(:, 2)] + ...
%          Az' * [g_t{n}(:, 1); g_t{n}(:, 2)];
% dgdtA = [dgdtA(1:GDense.J), dgdtA(GDense.J+1:end)];
Az_dense = FD_operator(G_dense, param.theta_z * (param.zmean - G_dense.z), param.sig_z*ones(G_dense.J,1), 2);

for n = 1:param.N
    
    % For KF, we are lucky and don't have to worry about BCs here because
    % drift is already inward-pointing.
    Aa_dense = FD_operator(G_dense, G.BH_dense * sim.s{n}, zeros(G_dense.J, 1), 1);

    AT = (Aa_dense + Az_dense)';
    
    if param.implicit_g
        % Implicit:
        B = 1/sim.dt(n) * speye(G_dense.J) - AT;
        b = sim.g{n} / sim.dt(n);
        sim.g{n+1} = B \ b;        
    else
        % Explicit: (for explicit, N should be >6 times larger than T)
        sim.g{n+1} = sim.g{n} + sim.dt(n) * AT * sim.g{n};
    end
    
    % Ensure positive density:
    sim.g{n+1}(sim.g{n+1} < 0) = 0;    
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - sum(sum(sim.g{n} * G_dense.dx ))) > 1e-6
        fprintf('KF not preserving mass.\n');
    end    
    
end


%% AGGREGATION & MARKET CLEARING
for n = 1:param.N    
    
    c_dense = G.BH_dense * sim.c{n};
    s_dense = G.BH_dense * sim.s{n};

    sim.mass(n,:) = sum(sim.g{n}*G_dense.dx);

    sim.B(n) = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.C(n) = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n) = sum(sum(s_dense .* sim.g{n} .* G_dense.dx));

end

sim.excess_bonds = sim.B;
sim.excess_goods = sim.Y - sim.C;
sim.excess_saving = sim.S;


%% COLLOCATION POINTS

DIFF_Y = interp1(sim.t, sim.excess_goods, param.nodes);
DIFF_B = interp1(sim.t, sim.excess_bonds, param.nodes);
DIFF_S = interp1(sim.t, sim.excess_saving, param.nodes);

diff = DIFF_S';


end




