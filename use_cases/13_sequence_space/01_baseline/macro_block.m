function sim = macro_block(x, z, ss, param)

%% SHOCK
sim.Z = param.Z * ones(param.N, 1);
sim.rho = param.rho * ones(param.N, 1);
sim.epsilon = param.epsilon * ones(param.N, 1);
sim.theta = zeros(param.N, 1);

switch param.shock_type
    case 'TFP'
        sim.Z = z;
        
    case 'demand'
        sim.rho = z;
        
    case 'cost-push'
        sim.epsilon = z;
    
    case 'monetary'
        sim.theta = z;
end


%% MACRO BLOCK

% Inputs:
sim.Y = reshape(x(1:param.N), [param.N, 1]);
sim.M = reshape(x(1+param.N:2*param.N), [param.N, 1]);

% Production and factor prices:
sim.N = sim.Y ./ sim.Z;
sim.w = sim.Z;

% Technology growth: (forward stencil)
sim.dZ = zeros(param.N, 1);
for n = 1:param.N-1
    sim.dZ(n) = (sim.Z(n+1) - sim.Z(n)) / sim.Z(n) / param.dt;
end
sim.dZ(param.N) = (ss.Z - sim.Z(param.N)) / sim.Z(param.N) / param.dt;

% Wage Phillips curve backward stencil: pi(n)-pi(n-1) = dt * X(n)
% sim.piw = zeros(param.N, 1);
% X = param.epsilon/param.delta * ((param.epsilon-1)/param.epsilon*(1+param.tau_L) * ss.w * ss.M ...
%     - param.v1(ss.N)) * ss.N;
% sim.piw(param.N) = ss.piw - param.dt * (param.rho * ss.piw + X);
% for n = param.N-1 : -1 : 1
%     X = sim.epsilon(n+1)/param.delta * (...
%           (sim.epsilon(n+1)-1)/sim.epsilon(n+1)*(1+param.tau_L) * sim.w(n+1) * sim.M(n+1) ...
%            - param.v1(sim.N(n+1)) ) * sim.N(n+1);
%     sim.piw(n) = sim.piw(n+1) - param.dt * (sim.rho(n+1) * sim.piw(n+1) + X);
% end

% Wage Phillips curve forward stencil: pi(n+1)-pi(n) = dt * X(n)
sim.piw = zeros(param.N, 1);
X = sim.epsilon(param.N)/param.delta * ((sim.epsilon(param.N)-1)/sim.epsilon(param.N)*(1+param.tau_L) ...
    * sim.w(param.N) * sim.M(param.N) - param.v1(sim.N(param.N))) * sim.N(param.N);
sim.piw(param.N) = (0 - param.dt * X) / (1 + param.dt * sim.rho(param.N));
for n = param.N-1 : -1 : 1
    X = sim.epsilon(n)/param.delta * ((sim.epsilon(n)-1)/sim.epsilon(n)*(1+param.tau_L) * sim.w(n) * sim.M(n) ...
        - param.v1(sim.N(n))) * sim.N(n);
    sim.piw(n) = (sim.piw(n+1) - param.dt * X) / (1 + param.dt * sim.rho(n));
end

% CPI inflation:
sim.pi = sim.piw - sim.dZ;

% Taylor rule:
sim.i = ss.r + param.lambda_pi * sim.pi + param.lambda_y * (sim.Y-ss.Y)/ss.Y + sim.theta;

% Fisher relation:
sim.r = sim.i - sim.pi;

end

