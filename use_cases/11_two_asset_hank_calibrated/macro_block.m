function sim = macro_block(x, z, ss, param)


%% SHOCK
sim.Z = param.Z * ones(param.N, 1);
sim.rho = param.rho * ones(param.N, 1);
% sim.epsilon = param.epsilon * ones(param.N, 1);
sim.G = param.G * ones(param.N, 1);
sim.eps = zeros(param.N, 1);

switch param.shock_type
    case 'monetary'
        sim.eps = z;
        
    case 'TFP'
        sim.Z = z;
        
    case 'demand'
        sim.rho = z;
        
    % case 'cost-push'
    %     sim.epsilon = z;
        
    case 'fiscal'
        sim.G  = z;
    
end


%% INPUTS
X = basis_fun_irf([], reshape(x, [1, numel(x)]), param.H, 3, param.bfun_type, param.t, "get_function");

sim.i  = X(:, 1);
sim.K  = X(:, 2);
sim.L  = X(:, 3);
% sim.Y = reshape(x(1+0*param.N : 1*param.N), [param.N, 1]);
% sim.I = reshape(x(1+1*param.N : 2*param.N), [param.N, 1]);
% sim.M = reshape(x(1+2*param.N : 3*param.N), [param.N, 1]);


%% MACRO BLOCK

% Production:
sim.Y  = sim.Z .* sim.K.^param.alpha .* sim.L.^(1-param.alpha);

% Taylor rule:
sim.pi = (sim.i - sim.eps - ss.r - param.lambda_Y * (sim.Y-ss.Y)/ss.Y) / param.lambda_pi;

% Fisher relation:
sim.r  = sim.i - sim.pi;

% NKPC: (no need to guess ik when assuming i is discount rate of firms)
% sim.P   = ones(param.N, 1);
% sim.dY  = zeros(param.N, 1);
% sim.dpi = zeros(param.N, 1);
% for n = 1:param.N-1
%     sim.P(n+1) = sim.P(n) .* (1+param.dt(n)*sim.pi(n));
%     sim.dY(n)  = (sim.Y(n+1) - sim.Y(n)) / param.dt(n);
%     sim.dpi(n) = (sim.pi(n+1) - sim.pi(n)) / param.dt(n);
% end
sim.mc = zeros(param.N, 1);
for n = 1:param.N-1
    sim.mc(n) = (param.epsilonF-1)/param.epsilonF + param.chiF/param.epsilonF ...
                 * (param.rho * sim.pi(n) - (sim.pi(n+1) - sim.pi(n)) / param.dt(n)) * sim.Y(n);
end
sim.mc(param.N) = (param.epsilonF-1)/param.epsilonF; % encodes terminal condition: pi(N) = 0

% Factor prices:
sim.rk = param.alpha * sim.mc .* sim.Y ./ sim.K;
% sim.ik = sim.rk .* sim.P;
sim.w  = 1/(1-param.tau_empl) * ( sim.Z .* sim.mc * (param.alpha^param.alpha * (1-param.alpha)^(1-param.alpha)) ...
         .* (sim.rk).^(-param.alpha) ).^(1./(1-param.alpha));

% Firm profits:
sim.Pi = (1-sim.mc) .* sim.Y;
% sim.L  = sim.Y .* (param.alpha/(1-param.alpha) .* (1-param.tau_empl) .* sim.w ./ sim.rk).^(-param.alpha);
% sim.K  = param.alpha./(1-param.alpha) .* (1-param.tau_empl) .* sim.L .* sim.w ./ sim.rk;
% sim.L = (1-param.alpha) * sim.mc .* sim.Y ./ ((1-param.tau_empl)*sim.w);
% sim.K = param.alpha * sim.mc .* sim.Y ./ sim.rk;

% Capital accumulation:
sim.dK = zeros(param.N, 1);
for n = 1:param.N-1
    sim.dK(n)  = (sim.K(n+1) - sim.K(n)) / param.dt(n);
end

% sim.K = ss.K * ones(param.N, 1);
% for n = 1:param.N-1
%     sim.K(n+1) = sim.K(n) + param.dt * (sim.I(n) - param.delta * sim.K(n));
% end
% % sim.I = sim.dK + param.delta * sim.K;
% sim.dK = sim.I - param.delta * sim.K;

% Investment and capital producer:
sim.gross_total_capital_accumulation = sim.dK + param.delta * sim.K;
sim.Q = param.solve_for_Q_from_cap_accumulation(sim.gross_total_capital_accumulation, sim.K);

sim.gross_total_investment_expenditure = param.gross_total_investment_expenditure(sim.Q, sim.K, 0);
sim.iotaQ = param.solve_for_iota(sim.Q);
assert( max(abs( sim.gross_total_capital_accumulation - param.gross_total_capital_accumulation(sim.Q, sim.K, 0) )) < 1e-8 );
% assert( max(abs( sim.Q - param.solveForQ(sim.iotaQ) )) < 1e-8 );
assert( max(abs( sim.iotaQ - sim.gross_total_capital_accumulation ./ sim.K )) < 1e-8 );

sim.I   = sim.gross_total_investment_expenditure;
sim.PiQ = param.PiQ(sim.Q, sim.K, zeros(param.N, 1));

% Hours and unemployment:
sim.U = param.la2 / (param.la1 + param.la2) * ones(param.N, 1);
sim.H = sim.L ./ (1-sim.U);

% Transfer:
sim.tau = sim.Pi + param.tau_lab * sim.w .* sim.L - param.UI * sim.U - sim.G - sim.r * param.gov_bond_supply;


end



