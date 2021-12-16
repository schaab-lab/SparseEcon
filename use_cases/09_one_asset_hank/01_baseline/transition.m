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
sim.l  = cell(param.N, 1);
sim.s  = cell(param.N, 1);
sim.sF = cell(param.N, 1);
sim.sB = cell(param.N, 1);

sim.NH = zeros(param.N, 1);
sim.C  = zeros(param.N, 1);
sim.S  = zeros(param.N, 1);


%% AGGREGATE TRANSITION PATH
X = basis_fun_irf([], reshape(PHI, [1, numel(PHI)]), param.H(1), param.H(2), param.bfun_type, sim.t, "get_function");

sim.N = X(:, 1); 
sim.i = X(:, 2);

sim.Z   = zeros(param.N, 1);
sim.eps = zeros(param.N, 1);
sim.G   = zeros(param.N, 1);
sim.rho = param.rho * ones(param.N, 1);
switch param.shock_type
    case 'monetary'
        sim.eps = shock;
        
    case 'productivity'
        sim.Z  = shock;
        
    case 'fiscal'
        sim.G  = shock;
        
    case 'demand'
        sim.rho = sim.rho + shock;
end

sim.Y = exp(sim.Z) .* sim.N;

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

sim.mc = (param.epsilon-1)/param.epsilon * (1 + param.chi/(param.epsilon-1) ...
         .* (sim.pi.*(sim.i - sim.pi - sim.dY./sim.Y) - sim.dpi));
sim.w  = 1/(1-param.tau_empl) * exp(sim.Z) .* sim.mc;

sim.Pi  = (1-sim.mc) .* sim.Y;
sim.U   = param.la2 / (param.la1 + param.la2);
sim.tau = sim.Pi + param.tau_lab * sim.w .* sim.N - sim.G - sim.r * param.gov_bond_supply;


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
    
    G.r = sim.r(n); G.w = sim.w(n); G.tau = sim.tau(n);
    
    % POLICY FUNCTIONS
    num0 = 1e-8; % numerical 0 for upwind scheme

    VaF = zeros(G.J, param.discrete_types);
    VaB = zeros(G.J, param.discrete_types);
    for j = 1:param.discrete_types
        VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    
    if n == param.N, c_r0 = G.r*G.a + (1-param.tau_lab)*G.w*param.zz + G.tau; c_l0 = c_r0; c0 = c_r0; end
    if n <  param.N, c_r0 = c_right; c_l0 = c_left; end
    
    f = @(x, a) x - G.tau - G.r*a - ((1-param.tau_lab)*G.w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta);
    f_fprime = @(x) 1 + param.gamma/param.eta*((1-param.tau_lab)*G.w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta-1);

    [c_right, ~] = newton_nonlin(f, f_fprime, c_r0, param.amax, param.crit);
    [c_left, ~] = newton_nonlin(f, f_fprime, c_l0, param.amin, param.crit);
    [c0, ~]      = newton_nonlin(f, f_fprime, c0,          G.a, param.crit);

    VaF(G.grid(:, 1)==1, :) = param.u1(c_right(G.grid(:, 1)==1, :));
    VaB(G.grid(:, 1)==0, :) = param.u1(c_left( G.grid(:, 1)==0, :));

    cF = param.u1inv(VaF);
    cB = param.u1inv(VaB);

    lF = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(cF));
    lB = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(cB));
    l0 = param.v1inv((1-param.tau_lab) * G.w * param.zz .* param.u1(c0));

    sF = G.r.*G.a + (1-param.tau_lab) * G.w.*param.zz.*lF + G.tau - cF;
    sB = G.r.*G.a + (1-param.tau_lab) * G.w.*param.zz.*lB + G.tau - cB;
    s0 = G.r.*G.a + (1-param.tau_lab) * G.w.*param.zz.*l0 + G.tau - c0;

    IF = (sF > num0);        % BC takes care of this: (G.grid(:, 1)<1)
    IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:, 1)>0)
    I0 = ~IF & ~IB;

    s = sF.*IF + sB.*IB;
    c = cF.*IF + cB.*IB + c0.*I0;
    l = lF.*IF + lB.*IB + l0.*I0;
    u = param.u(c) - param.v(l);

    % CONSTRUCT FD OPERATORS
    Aa{1} = FD_operator(G, s(:, 1), zeros(G.J, 1), 1, '1');
    Aa{2} = FD_operator(G, s(:, 2), zeros(G.J, 1), 1, '2');

    A = blkdiag(Aa{1}, Aa{2}) + Az;

    B = (1/sim.dt(n) + sim.rho(n))*speye(2*G.J) - A;
    b = [u(:, 1); u(:, 2)] + [V(:, 1); V(:, 2)] / sim.dt(n);

    % SOLVE LINEAR SYSTEM
    V_new = B\b;
    V = [V_new(1:G.J), V_new(1+G.J:end)];
    
    if ~isreal(V), disp('Complex values detected!'); diff = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n}=V; sim.u{n}=u; sim.c{n}=c; sim.s{n}=s; sim.l{n}=l;
end


%% SOLVE KF FORWARDS
Az_dense = [-speye(G_dense.J)*param.la1,  speye(G_dense.J)*param.la1; ...
             speye(G_dense.J)*param.la2, -speye(G_dense.J)*param.la2];

for n = 1:param.N
    
    % For KF, we are lucky and don't have to worry about BCs here because
    % drift is already inward-pointing. 
    Aa_dense{1} = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 1), zeros(G_dense.J, 1), 1, '1');
    Aa_dense{2} = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 2), zeros(G_dense.J, 1), 1, '2');

    AT = (blkdiag(Aa_dense{1}, Aa_dense{2}) + Az_dense)';
    
    if param.implicit_g
        % Implicit:
        B = 1/sim.dt(n) * speye(param.discrete_types*G_dense.J) - AT;
        b = [sim.g{n}(:, 1); sim.g{n}(:, 2)]/sim.dt(n);

        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    elseif ~param.implicit_g && param.reso_sim_KF==1
        % Explicit: (for explicit, N must be >6 times larger than T)
        gg = [sim.g{n}(:, 1); sim.g{n}(:, 2)] + sim.dt(n) * AT * [sim.g{n}(:, 1); sim.g{n}(:, 2)];
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    elseif ~param.implicit_g && param.reso_sim_KF>1
        if n<param.N
            t_KF  = linspace(sim.t(n), sim.t(n+1), param.reso_sim_KF+1);
            dt_KF = t_KF(2) - t_KF(1);
        end
        g_KF    = cell(1, param.reso_sim_KF+1);
        g_KF{1} = sim.g{n};
        
        for m = 1:param.reso_sim_KF
            dgdtA = [Aa_dense{1}' * g_KF{m}(:, 1); Aa_dense{2}' * g_KF{m}(:, 2)] + ...                    
                     Az_dense' * [g_KF{m}(:, 1); g_KF{m}(:, 2)];
            dgdtA = [dgdtA(1:G_dense.J), dgdtA(G_dense.J+1:end)];

            g_KF{m+1} = g_KF{m} + dgdtA * dt_KF;
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


%% AGGREGATION & MARKET CLEARING
for n = 1:param.N    
    
    l_dense = G.BH_dense * sim.l{n};
    c_dense = G.BH_dense * sim.c{n};
    s_dense = G.BH_dense * sim.s{n};

    sim.mass(n, :) = sum(sim.g{n}*G_dense.dx);

    sim.LH(n) = sum(sum(l_dense .* param.zz .* sim.g{n} .* G_dense.dx));
    sim.BH(n) = sum(sum(G_dense.a .* sim.g{n} .* G_dense.dx));
    sim.C(n)  = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n)  = sum(sum(s_dense .* sim.g{n} .* G_dense.dx));

end

sim.excess_bonds  = sim.BH - param.gov_bond_supply;
sim.excess_saving = sim.S;
sim.excess_labor  = sim.NH - sim.N;
sim.excess_goods  = sim.Y - sim.C;


%% COLLOCATION POINTS

DIFF_Y = interp1(sim.t, sim.excess_goods,  param.nodes);
DIFF_L = interp1(sim.t, sim.excess_labor,  param.nodes);
DIFF_B = interp1(sim.t, sim.excess_bonds,  param.nodes);
DIFF_S = interp1(sim.t, sim.excess_saving, param.nodes);

diff = [DIFF_Y, DIFF_S]';


end




