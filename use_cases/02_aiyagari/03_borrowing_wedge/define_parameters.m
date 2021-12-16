function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 2; param.surplus = [3, 0];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [7, 2]; % vector of "surplus" for dense grid

param.amin = -2;
param.amax = 50;
param.zmin = 0.3;
param.zmax = 1.5;

param.min = [param.amin, param.zmin];
param.max = [param.amax, param.zmax];

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 20;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_tol\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 100;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho   = 0.05;
param.gamma = 2;
param.alpha = 0.33;
param.delta = 0.05;
param.wedge = 0.01;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

param.v     = @(x) x.^(1+param.eta) / (1+param.eta);
param.v1    = @(x) x.^param.eta;
param.v1inv = @(x) x.^(1/param.eta);

% Earnings parameters:
param.zmean   = (param.zmax + param.zmin)/2;
param.theta_z = 0.25;
param.sig_z   = 0.02;
param.L       = param.zmean;

% Government parameters:
param.gov_bond_supply = 0;


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