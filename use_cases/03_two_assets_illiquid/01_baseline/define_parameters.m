function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
% param.l = 3; param.surplus = [3,3];
param.l = 0; param.surplus = [4,4];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [4,4]; % vector of "surplus" for dense grid

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
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_told\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 300;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 1000;
param.crit_KF  = 1e-7;


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho       = 0.035;%0.012728925130000;%0.05;
param.deathrate = 0.005555555555556;
param.gamma     = 2;
param.eta       = 2;
param.delta     = 0.0125;

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

% Earnings parameters:
param.z1 = 0;
param.z2 = 1;

param.U0 = 0.089;

param.la1 = 0.97; % 1/x implies x quarters of duration
param.la2 = 0.09; %param.U0 * param.la1 / (1- param.U0);

% param.L  = param.la1/(param.la1+param.la2) * exp(param.z1) + param.la2/(param.la1+param.la2) * exp(param.z2);
param.L  = param.la2/(param.la1+param.la2) * param.z1 + param.la1/(param.la1+param.la2) * param.z2;
param.zz = [param.z1, param.z2];
param.discrete_types = 2;

% Portfolio choice parameters:
param.dmax = 1e+10;
param.psi0 = 0.04383; %0;
param.psi1 = 1 - 0.04383; %0.5;
param.psi2 = 1 + 0.40176; %2;
param.psi3 = 0.03*2.92/4;

param.cap_adj_model = 'KMV_delta_offset';
param.ZQmean = 0;
param.xi = 1;

% Firms and Inflation:
param.alpha = 0.38;
param.epsilonF = 10;
param.chiF  = 100;
param.kappa = 100;

% Unions:
param.wage_rigidity_estimation = 'ACEL';
param.epsilonW = 21;
RHS = @(theta,rho) (1-theta) * (1-theta * 1/(1+rho)) / theta;
chi = @(epsilonW,theta,rho) (epsilonW-1) / RHS(theta, rho);
switch param.wage_rigidity_estimation
    case 'CEE'
        theta = 0.64;
    case 'ACEL' 
        theta = 0.78;
end
param.chi = chi(param.epsilonW, theta, param.rho);

% Government:
param.lambda_pi = 1.05;
param.lambda_Y  = 0.02;

param.tau_lump = 0;
param.tau_lab = 0.20;
param.UI = 0.20;
param.G = 0;
param.gov_bond_supply = 1.00;

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


        param.gross_total_capital_accumulation   = @(Q, K, ZQ) param.Phi(param.cap_investment(Q,ZQ), ZQ);
        param.gross_total_investment_expenditure = @(Q, K, ZQ) param.cap_investment(Q,ZQ);
        param.PiQ = @(Q, K, ZQ) Q .* param.Phi(param.cap_investment(Q,ZQ),ZQ) - param.cap_investment(Q,ZQ);
        
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
        param.PiQ = @(Q, K, ZQ) Q.*param.solve_for_iota(Q).*K - param.solve_for_iota(Q).*K ...
                                                    - param.Phi(param.solve_for_iota(Q)).*K;
end


%% VARIABLE INPUTS

% Parse inputs:
p = inputParser;
p.CaseSensitive = true;
for f = fieldnames(param)'
    p.addParameter(f{:}, param.(f{:}));
end
parse(p, varargin{:});
param = p.Results;

% Update parameters


end

