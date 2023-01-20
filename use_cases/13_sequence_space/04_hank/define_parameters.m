function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 6;
param.d = 1; param.d_idio = 1; param.d_agg = 0;

param.l_dense = 6; % vector of "surplus" for dense grid

param.amin = -1;
param.amax = 15;

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 1;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_tol\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 100;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;


%% JACOBIAN TUNING PARAMETERS
param.phi_jacobian = 1000;
param.psi_jacobian = 0;


%% TRANSITION DYNAMICS PARAMETERS
param.T = 100; 
param.N = 120;

param.implicit_g = 0;


%% ECONOMIC PARAMETERS

% Shock:
param.shock_type = 'TFP';
param.shock_percent = 0.01;
param.shock_theta = log(2);

% Household parameters:
param.rho = 0.02;
param.eta = 2;
param.gamma = 2;
param.delta = 100;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);
param.u2    = @(x) -param.gamma * x.^(-param.gamma-1);

param.v     = @(x) x.^(1+param.eta) / (1+param.eta); 
param.v1    = @(x) x.^param.eta;
param.v1inv = @(x) x.^(1/param.eta);

% Earnings parameters:
param.zz  = [0.8, 1.2];
param.la1 = 1/3;
param.la2 = 1/3;
param.L   = param.la2/(param.la1+param.la2) * param.zz(1) + param.la1/(param.la1+param.la2) * param.zz(2);

param.discrete_types = numel(param.zz);

% TFP:
param.Z = 1;

% Unions:
param.epsilon = 10;

% Government:
param.tau_L = 0;
param.lambda_pi = 1.5;
param.lambda_y = 0;


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
param.min = param.amin;
param.max = param.amax;

% Transition path:
param.t = linspace(0, param.T, param.N)';
param.dt = param.t(2) - param.t(1); %diff(param.t); param.dt(param.N) = param.dt(param.N-1);

% Earnings: 
param.L   = param.la2/(param.la1+param.la2) * param.zz(1) + param.la1/(param.la1+param.la2) * param.zz(2);

% Unions:
param.kappa = param.epsilon/param.delta * (param.epsilon-1)/param.epsilon*(1+param.tau_L);

% Shocks:
switch param.shock_type
    case 'TFP'
        param.shock_level = param.shock_percent * param.Z;
        % param.shock_theta = log(2);

    case 'demand'
        param.shock_level = 0.25 * param.rho;
        % param.shock_theta = log(2);
        
    case 'cost-push'
        param.shock_level = 0.10 * param.epsilon;
        % param.shock_theta = log(2);
        
end

end
