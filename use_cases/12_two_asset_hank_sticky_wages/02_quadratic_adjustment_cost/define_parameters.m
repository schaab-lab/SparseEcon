function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
% param.l = 3; param.surplus = [3, 3];
param.l = 0; param.surplus = [5, 5];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [4, 4]; % vector of "surplus" for dense grid

param.amin = -1;
param.amax = 7;
param.kmin = 0; %0.001;
param.kmax = 50;

param.min = [param.amin, param.kmin];
param.max = [param.amax, param.kmax];

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 20;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_tol\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 300;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-7;


%% JACOBIAN TUNING PARAMETERS
param.phi_jacobian = 1000;
param.psi_jacobian = 0;


%% TRANSITION DYNAMICS PARAMETERS
param.T = 100; 
param.N = 120;

param.implicit_g = 0;


%% ECONOMIC PARAMETERS

% Model:
param.natural = 0;
param.compute_natural = 0;
% param.compute_no_policy = 0;
% param.production_run = 0;
% param.discretion_no_inflation_penalty = 0;
% param.discretion_no_distributional_penalty = 0;

% Shock:
param.shock_type = 'TFP';
param.shock_percent = 0.01;
param.shock_theta = log(2);

% Household parameters:
param.rho   = 0.035; %0.075;%0.012728925130000;%0.05;
param.gamma = 2;
param.eta   = 2;
param.delta = 0.03;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);
param.u2    = @(x) -param.gamma * x.^(-param.gamma-1);
param.u3    = @(x) param.gamma * (1+param.gamma) * x.^(-param.gamma-2);

param.v     = @(x) x.^(1+param.eta) / (1+param.eta); 
param.v1    = @(x) x.^param.eta;
param.v1inv = @(x) x.^(1/param.eta);
param.v2    = @(x) param.eta * x.^(param.eta-1);

% Earnings parameters:
param.zz  = [0.8, 1.2];
param.la1 = 1/3;
param.la2 = 1/3;
param.L   = param.la2/(param.la1+param.la2) * param.zz(1) + param.la1/(param.la1+param.la2) * param.zz(2);

% param.z1 = 0;
% param.z2 = 1;
% param.U0 = 0.089;

% param.la1 = 0.97; % 1/x implies x quarters of duration
% param.la2 = 0.09; %param.U0 * param.la1 / (1- param.U0);
% param.L  = param.la2/(param.la1+param.la2) * param.z1 + param.la1/(param.la1+param.la2) * param.z2;
% param.zz = [param.z1, param.z2];
param.discrete_types = 2;

% Household portfolio choice:
param.dmax = 1e+10;
param.psi0 = 0.04383; %0;
param.psi1 = 1 - 0.04383; %0.5;
param.psi2 = 1 + 0.40176; %2;
param.psi3 = 0.03*2.92/4;

param.cap_adj_model = 'KMV_delta_offset';
param.ZQmean = 0;
param.xi = 1;

% TFP:
param.Z = 1;

% Firms and Inflation:
param.alpha = 0.38;
param.kappa = 100;

% Unions:
param.epsilon = 10;
param.chi = 100;

% Government:
param.policy_shock = 0.002;
param.theta_policy = log(2);

param.lambda_pi = 1.5;
param.lambda_y = 0;

param.tau_L = 0;
param.tau_lab = 0.20;

param.G = 0;
param.gov_bond_supply = 1.00;


%% VARIABLE INPUTS

% Parse inputs:
p = inputParser;
p.CaseSensitive = true;
for f = fieldnames(param)'
    p.addParameter(f{:}, param.(f{:}));
end
parse(p, varargin{:});
param = p.Results;

%% UPDATE PARAMETERS

% Grid:
param.min = [param.amin, param.kmin];
param.max = [param.amax, param.kmax];

% Transition path:
param.t = linspace(0, param.T, param.N)';
param.dt = param.t(2) - param.t(1); %diff(param.t); param.dt(param.N) = param.dt(param.N-1);

% Household preferences:
param.v  = @(h) h.^(1+param.eta)/(1+param.eta);
param.v1 = @(h) h.^param.eta;

if param.gamma == 1
    param.u     = @(c) log(c);
    param.uinv  = @(u) exp(u);
    param.u1    = @(c) 1./c;
    param.u1inv = @(u) 1./u;
    param.u2    = @(c) - 1 ./ (c.^2);
else 
    param.u     = @(c) c.^(1-param.gamma) ./ (1-param.gamma);
    param.uinv  = @(u) (1-param.gamma) * u.^( 1 / (1-param.gamma));
    param.u1    = @(c) c.^(-param.gamma);
    param.u1inv = @(u) u.^(-1/param.gamma);
    param.u2    = @(c) -param.gamma * c.^(-param.gamma - 1); 
end

% Earnings: 
param.L   = param.la2/(param.la1+param.la2) * param.zz(1) + param.la1/(param.la1+param.la2) * param.zz(2);

% Aggregate investment adjustment cost:
switch param.cap_adj_model
    % BRU-SAN model of adjustment cost:
    case 'BS'
        param.Phi          = @(iota, ZQ) exp(ZQ) ./ param.kappa .* (sqrt(1 + 2*param.kappa*iota) - 1);
        param.Phi_prime    = @(iota, ZQ) exp(ZQ) ./ ( sqrt( 1 + 2*param.kappa * iota ));
        param.Phi_inv      = @(iota, ZQ) ( (1+param.kappa*iota).^2 - 1 ) / (2 * param.kappa);

        param.cap_investment              = @(Q, ZQ) ( (exp(ZQ).*Q) .^2 - 1) ./ (2 * param.kappa);
        param.cap_investment_prime        = @(Q, ZQ) (exp(ZQ)).^2 / param.kappa .* Q;
        param.cap_investment_double_prime = @(Q, ZQ) (exp(ZQ)).^2 / param.kappa;


        param.gross_total_capital_accumulation   = @(Q, K, ZQ) param.Phi(param.cap_investment(Q, ZQ), ZQ);
        param.gross_total_investment_expenditure = @(Q, K, ZQ) param.cap_investment(Q, ZQ);
        param.PiQ = @(Q, K, ZQ) Q .* param.Phi(param.cap_investment(Q, ZQ), ZQ) - param.cap_investment(Q, ZQ);
        
    % KMV model of adjustment cost:
    case 'KMV'
        param.Phi = @(iota) param.kappa / 2 * (iota).^2;
        param.Phi_prime = @(iota) param.kappa * iota;

        param.solve_for_Q = @(iota) 1 + param.Phi_prime(iota);
        param.solve_for_iota = @(Q) (Q-1) / param.kappa; 

        param.gross_total_capital_accumulation   = @(Q, K, ZQ) param.solve_for_iota(Q) .* K;
        param.gross_total_investment_expenditure = @(Q, K, ZQ) param.solve_for_iota(Q) .* K ...
                                                    + param.Phi(param.solve_for_iota(Q)).*K;
        param.PiQ = @(Q, K, ZQ) Q.*param.solve_for_iota(Q).*K - param.solve_for_iota(Q).*K ...
                                                    - param.Phi(param.solve_for_iota(Q)).*K;

    case 'KMV_delta_offset'
        param.Phi = @(iota) param.kappa / 2 * (iota - param.delta).^2;
        param.Phi_prime = @(iota) param.kappa * (iota - param.delta);

        param.solve_for_Q = @(iota) 1 + param.Phi_prime(iota);
        param.solve_for_iota = @(Q) (Q-1) / param.kappa + param.delta; 

        param.gross_total_capital_accumulation   = @(Q, K, ZQ) param.solve_for_iota(Q) .* K;
        param.gross_total_investment_expenditure = @(Q, K, ZQ) param.solve_for_iota(Q) .* K ...
                                                    + param.Phi(param.solve_for_iota(Q)).*K;

        param.solve_for_Q_from_cap_accumulation = @(X, K) 1 + param.kappa * (X ./ K - param.delta);
        
        param.PiQ = @(Q, K, ZQ) Q.*param.solve_for_iota(Q).*K - param.solve_for_iota(Q).*K ...
                                                    - param.Phi(param.solve_for_iota(Q)).*K;
end

% Shocks:
switch param.shock_type
    case 'TFP'
        param.shock_level = param.shock_percent * param.Z;
        % param.shock_theta = log(2);

    case 'demand'
        param.shock_level = 0.05 * param.rho;
        % param.shock_level = 0.25 * param.rho;
        % param.shock_theta = log(2);
        
    case 'cost-push'
        param.shock_level = 0.10 * param.epsilon;
        % param.shock_theta = log(2);
        
end

end

