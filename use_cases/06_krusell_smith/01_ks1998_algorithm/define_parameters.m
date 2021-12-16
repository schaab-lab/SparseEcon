function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 2;
param.d = 3; param.d_idio = 1; param.d_agg = 2; % (a,Z,K)
param.surplus = [5, 0, 2];

param.l_dense = 6; % vector of "surplus" for dense grid

param.amin = 0;
param.amax = 20;

param.Zmin = -0.10;
param.Zmax =  0.10;

param.min = [param.amin, param.Zmin];
param.max = [param.amax, param.Zmax];

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 20;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_told\n'); end

% KS algorithm:
param.estimation_model_type = 'linear';
param.lambda_LOM = 0.80;
param.max_KS = 100;
param.crit_KS = 1e-6;

% Simulation parameters:
param.T  = 200;
param.N  = 2000;
param.t  = linspace(0, param.T, param.N);
param.dt = param.t(2) - param.t(1);

param.burn_in = 50;
param.n_data  = (round(param.burn_in/param.dt):param.N)';

% PDE tuning parameters:
param.Delta = 1000;
param.maxit = 300;
param.crit  = 1e-7;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho = 0.05;
param.gamma = 2;
param.alpha = 1/3;
param.delta = 0.05;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

% Earnings parameters:
param.zz  = [0.9, 1.1];
param.la1 = 1/3;
param.la2 = 1/3;
param.L   = param.la1/(param.la1+param.la2) * param.zz(1) + param.la2/(param.la1+param.la2) * param.zz(2);

param.discrete_types = numel(param.zz);

% Aggregate risk:
param.Zmean = 0;
param.thetaZ = 0.25;
param.sigmaZ = 0.007;


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
param.t  = linspace(0, param.T, param.N);
param.dt = param.t(2) - param.t(1);

param.n_data = (round(param.burn_in/param.dt):param.N)';

end