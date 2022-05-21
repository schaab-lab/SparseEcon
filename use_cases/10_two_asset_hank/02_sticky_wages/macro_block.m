function sim = macro_block(x, t, z, ss, param)

%% SHOCK
sim.Z = param.Z * ones(param.N, 1);
sim.rho = param.rho * ones(param.N, 1);
sim.epsilon = param.epsilon * ones(param.N, 1);
sim.G = param.G * ones(param.N, 1);

switch param.shock_type
    case 'TFP'
        sim.Z = z;
        
    case 'demand'
        sim.rho = z;
        
    case 'cost-push'
        sim.epsilon = z;
        
end


%% NATURAL ALLOCATION
% if param.natural
%     sim.r = reshape(x(1:param.N), [param.N, 1]);
%     sim.M = reshape(x(1+param.N:2*param.N), [param.N, 1]);
% 
%     sim.w = sim.Z;
%     
%     sim.N = param.v1inv((sim.epsilon-1)./sim.epsilon*(1+param.tau_L) .* sim.w .* sim.M);
%     sim.Y = sim.Z .* sim.N;
%     
%     sim.piw = zeros(param.N, 1);
% end


%% STANDARD MACRO BLOCK
% if ~param.natural

% Inputs:
sim.Y = reshape(x(1+0*param.N : 1*param.N), [param.N, 1]);
sim.K = reshape(x(1+1*param.N : 2*param.N), [param.N, 1]);
sim.M = reshape(x(1+2*param.N : 3*param.N), [param.N, 1]);
sim.theta = t;

% Production and factor prices:
sim.N  = (sim.Y ./ (sim.Z .* sim.K.^param.alpha)) .^ (1/(1-param.alpha));
sim.rk = param.alpha * sim.Y ./ sim.K;
sim.w  = (1-param.alpha) * sim.Y ./ sim.N;

% Growth: (forward stencil)
sim.dZ = zeros(param.N, 1);
sim.dY = zeros(param.N, 1);
sim.dN = zeros(param.N, 1);
for n = 1:param.N-1
    sim.dZ(n) = (sim.Z(n+1) - sim.Z(n)) / sim.Z(n) / param.dt;
    sim.dY(n) = (sim.Y(n+1) - sim.Y(n)) / sim.Y(n) / param.dt;
    sim.dN(n) = (sim.N(n+1) - sim.N(n)) / sim.N(n) / param.dt;
end
sim.dZ(param.N) = (ss.Z - sim.Z(param.N)) / sim.Z(param.N) / param.dt;
sim.dY(param.N) = (ss.Y - sim.Y(param.N)) / sim.Y(param.N) / param.dt;
sim.dN(param.N) = (ss.N - sim.N(param.N)) / sim.N(param.N) / param.dt;

% Wage Phillips curve backward stencil: pi(n)-pi(n-1) = dt * X(n)
% sim.piw = zeros(param.N, 1);
% X = param.epsilon/param.chi * ((param.epsilon-1)/param.epsilon*(1+param.tau_L) * ss.w * ss.M ...
%     - param.v1(ss.N)) * ss.N;
% sim.piw(param.N) = ss.piw - param.dt * (param.rho * ss.piw + X);
% for n = param.N-1 : -1 : 1
%     X = sim.epsilon(n+1)/param.chi * (...
%           (sim.epsilon(n+1)-1)/sim.epsilon(n+1)*(1+param.tau_L) * sim.w(n+1) * sim.M(n+1) ...
%            - param.v1(sim.N(n+1)) ) * sim.N(n+1);
%     sim.piw(n) = sim.piw(n+1) - param.dt * (sim.rho(n+1) * sim.piw(n+1) + X);
% end

% Wage Phillips curve forward stencil: pi(n+1)-pi(n) = dt * X(n)
sim.piw = zeros(param.N, 1);
X = sim.epsilon(param.N)/param.chi * ((sim.epsilon(param.N)-1)/sim.epsilon(param.N)*(1+param.tau_L)*(1-param.tau_lab) ...
    * sim.w(param.N) * sim.M(param.N) - param.v1(sim.N(param.N))) * sim.N(param.N);
sim.piw(param.N) = (ss.piw - param.dt * X) / (1 + param.dt * sim.rho(param.N));
for n = param.N-1 : -1 : 1
    X = sim.epsilon(n)/param.chi * ((sim.epsilon(n)-1)/sim.epsilon(n)*(1+param.tau_L)*(1-param.tau_lab) * sim.w(n) * sim.M(n) ...
        - param.v1(sim.N(n))) * sim.N(n);
    sim.piw(n) = (sim.piw(n+1) - param.dt * X) / (1 + param.dt * sim.rho(n));
end

% Capital production:
sim.dK = zeros(param.N, 1);
for n = 1:param.N-1
sim.dK(n)  = (sim.K(n+1) - sim.K(n)) / param.dt;
end

sim.gross_total_capital_accumulation = sim.dK + param.delta * sim.K;
sim.Q = param.solve_for_Q_from_cap_accumulation(sim.gross_total_capital_accumulation, sim.K);

sim.gross_total_investment_expenditure = param.gross_total_investment_expenditure(sim.Q, sim.K, 0);
sim.iotaQ = param.solve_for_iota(sim.Q);
assert( max(abs( sim.gross_total_capital_accumulation - param.gross_total_capital_accumulation(sim.Q, sim.K, 0) )) < 1e-8 );
assert( max(abs( sim.Q - param.solve_for_Q(sim.iotaQ) )) < 1e-8 );
assert( max(abs( sim.iotaQ - sim.gross_total_capital_accumulation ./ sim.K )) < 1e-8 );

sim.I   = sim.gross_total_investment_expenditure;
sim.PiQ = param.PiQ(sim.Q, sim.K, zeros(param.N, 1));

% CPI inflation:
sim.pi = sim.piw - (sim.dY - sim.dN);

% Taylor rule:
sim.i = (ss.r + ss.piw) + param.lambda_pi * sim.pi + param.lambda_y * (sim.Y-ss.Y)/ss.Y + sim.theta;

% Fisher relation:
sim.r = sim.i - sim.pi;

% Transfer:
sim.tau = param.tau_lab * sim.w .* sim.N - sim.G - sim.r * param.gov_bond_supply;

% end
end



