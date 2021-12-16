function [diff, G, G_dense, sim] = transition(PHI, G, G_dense, shock, ss, param)

%% SETUP
sim.N  = param.N;
sim.t  = param.t;
sim.dt = param.dt;

% PREALLOCATE
sim.V  = cell(param.N,1);
sim.A  = cell(param.N,1);
sim.g  = cell(param.N,1);
sim.c  = cell(param.N,1);
sim.l  = cell(param.N,1);
sim.s  = cell(param.N,1);
sim.sF = cell(param.N,1);
sim.sB = cell(param.N,1);

sim.KH = zeros(param.N,1);
sim.LH = zeros(param.N,1);
sim.C  = zeros(param.N,1);
sim.S  = zeros(param.N,1);
sim.I  = zeros(param.N,1);


%% AGGREGATE TRANSITION PATH
X = basis_fun_irf([], reshape(PHI, [1, numel(PHI)]), param.H(1), param.H(2), ...
    param.bfun_type, sim.t, "get_function");

sim.Y = X(:,1); 
sim.L = X(:,2);
sim.Z = shock;

sim.K  = (sim.Y ./ (exp(sim.Z) .* sim.L.^(1-param.alpha)) ).^(1/param.alpha);
sim.w  = (1-param.alpha) * sim.Y ./ sim.L;
sim.rk = param.alpha * sim.Y ./ sim.K;
sim.r  = sim.rk - param.delta; 


%% SOLVE VFI BACKWARDS
V = ss.V;
sim.g{1} = ss.g;

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

for n = param.N:-1:1
    
    G.r = sim.r(n); G.w = sim.w(n);
    
    % POLICY FUNCTIONS
    num0 = 1e-8; % numerical 0 for upwind scheme

    VkF = zeros(G.J, param.discrete_types);
    VkB = zeros(G.J, param.discrete_types);
    for j = 1:param.discrete_types
        VkF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VkB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    
    if n == param.N, c_r0 = G.r * G.k + G.w * param.zz; c_l0 = c_r0; c0 = c_r0; end
    if n <  param.N, c_r0 = c_right; c_l0 = c_left; end
    
    f = @(x,k) x - G.r*k - (G.w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta);
    f_prime = @(x) 1 + param.gamma/param.eta*(G.w*param.zz).^(1/param.eta + 1) .* x.^(-param.gamma/param.eta-1);

    [c_right, ~] = newton_nonlin(f, f_prime, c_r0, param.kmax, param.crit);
    [c_left,  ~] = newton_nonlin(f, f_prime, c_l0, param.kmin, param.crit);
    [c0, ~]      = newton_nonlin(f, f_prime, c0,          G.k, param.crit);

    VkF(G.grid(:, 1) == 1, :) = param.u1(c_right(G.grid(:, 1) == 1, :));
    VkB(G.grid(:, 1) == 0, :) = param.u1(c_left( G.grid(:, 1) == 0, :));
    
    cF = param.u1inv(VkF);
    cB = param.u1inv(VkB);
    
    lF = param.v1inv(G.w * param.zz .* param.u1(cF));
    lB = param.v1inv(G.w * param.zz .* param.u1(cB));
    l0 = param.v1inv(G.w * param.zz .* param.u1(c0));
    
    sF = G.r.*G.k + G.w.*param.zz.*lF - cF;
    sB = G.r.*G.k + G.w.*param.zz.*lB - cB;
    s0 = G.r.*G.k + G.w.*param.zz.*l0 - c0;
    
    IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
    IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
    I0 = ~IF & ~IB;
    
    s = sF.*IF + sB.*IB;
    c = cF.*IF + cB.*IB + c0.*I0;
    l = lF.*IF + lB.*IB + l0.*I0;
    u = param.u(c) - param.v(l);
    
    % CONSTRUCT FD OPERATORS
    Ak{1} = FD_operator(G, s(:, 1), zeros(G.J, 1), 1, '1');
    Ak{2} = FD_operator(G, s(:, 2), zeros(G.J, 1), 1, '2');
    
    A = blkdiag(Ak{1}, Ak{2}) + Az;
    
    B = (1/sim.dt(n) + param.rho)*speye(2*G.J) - A;
    b = [u(:, 1); u(:, 2)] + [V(:, 1); V(:, 2)] / sim.dt(n);
    
    % SOLVE LINEAR SYSTEM
    V_new = B\b;
    V = [V_new(1:G.J), V_new(1+G.J:end)];
    
    if ~isreal(V), disp('Complex values detected!'); diff = NaN(1); return; end
    
    % RECORD DATA
    sim.V{n} = V; sim.u{n} = u; sim.c{n} = c; sim.s{n} = s; sim.l{n} = l;
end


%% SOLVE KF FORWARDS
% dgdtA = [Aa1' * g_t{n}(:, 1); Aa2' * g_t{n}(:, 2)] + ...
%         [Ak1' * g_t{n}(:, 1); Ak2' * g_t{n}(:, 2)] + ...
%          Az' * [g_t{n}(:, 1); g_t{n}(:, 2)];
% dgdtA = [dgdtA(1:GDense.J), dgdtA(GDense.J+1:end)];
Az_dense = [-speye(G_dense.J)*param.la1,  speye(G_dense.J)*param.la1; ...
             speye(G_dense.J)*param.la2, -speye(G_dense.J)*param.la2];

for n = 1:param.N
    
    % For KF, we are lucky and don't have to worry about BCs here because
    % drift is already inward-pointing. 
    Ak_dense{1} = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 1), zeros(G_dense.J, 1), 1, '1');
    Ak_dense{2} = FD_operator(G_dense, G.BH_dense * sim.s{n}(:, 2), zeros(G_dense.J, 1), 1, '2');

    AT = (blkdiag(Ak_dense{1}, Ak_dense{2}) + Az_dense)';
    
    if param.implicit_g
        % Implicit:
        B = 1/sim.dt(n) * speye(param.discrete_types*G_dense.J) - AT;
        b = [sim.g{n}(:, 1); sim.g{n}(:, 2)]/sim.dt(n);

        gg = B \ b;
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    else
        % Explicit: (for explicit, N must be >6 times larger than T)
        gg = [sim.g{n}(:, 1); sim.g{n}(:, 2)] + sim.dt(n) * AT * [sim.g{n}(:, 1); sim.g{n}(:, 2)];
        sim.g{n+1} = [gg(1:G_dense.J), gg(1+G_dense.J:end)];
    end
    
    for j = 1:param.discrete_types, sim.g{n+1}(sim.g{n+1}(:, j) < 0, j) = 0; end
    
    if abs( sum(sum(sim.g{n+1} * G_dense.dx )) - sum(sum(sim.g{n} * G_dense.dx ))) > 1e-6
        fprintf('KF not preserving mass.\n');
    end    
    
end


%% AGGREGATION & MARKET CLEARING
for n = 1:param.N    
    
    l_dense = G.BH_dense * sim.l{n};
    c_dense = G.BH_dense * sim.c{n};
    s_dense = G.BH_dense * sim.s{n};

    sim.mass(n,:) = sum(sim.g{n}*G_dense.dx);

    sim.LH(n) = sum(sum(l_dense .* param.zz .* sim.g{n} .* G_dense.dx));
    sim.KH(n) = sum(sum(G_dense.k .* sim.g{n} .* G_dense.dx));
    sim.C(n)  = sum(sum(c_dense .* sim.g{n} .* G_dense.dx));
    sim.S(n)  = sum(sum(s_dense .* sim.g{n} .* G_dense.dx));
    sim.I(n)  = sim.S(n) + param.delta * sim.KH(n); 

end

sim.excess_capital = sim.KH - sim.K;
sim.excess_labor   = sim.LH - sim.L;
sim.excess_goods   = sim.Y - sim.C - sim.I;


%% COLLOCATION POINTS

DIFF_Y = interp1(sim.t, sim.excess_goods,   param.nodes);
DIFF_L = interp1(sim.t, sim.excess_labor,   param.nodes);
DIFF_K = interp1(sim.t, sim.excess_capital, param.nodes);

diff = [DIFF_Y, DIFF_L]';


end




