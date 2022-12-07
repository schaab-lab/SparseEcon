function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 2; param.surplus = [2, 4];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [6, 3]; % vector of "surplus" for dense grid

param.amin = 0;
param.amax = 90;
param.zmin = 0.3;
param.zmax = 2.2;

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
param.rho   = 0.07;
param.gamma = 2;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

% Earnings parameters:
param.discrete_types = 1; %numel(param.zz);
param.L     = 1;
param.theta = 0.02;

param.zmean   = 0;
param.theta_z = 0.02;
param.sig_z   = 0.1265;

% Production parameters
param.delta = 0.05;

% fixed cost input capital
param.fkU = 0;
param.fkP = 10;

% fixed cost input labor
param.flU = 0;
param.flP = 8;

% fixed cost to operate technology
param.fyU = 0;
param.fyP = 0;

% DRS
param.RS    = 0.75;
param.eta   = 0.4;
param.alpha = param.eta * param.RS;
param.beta  = (1-param.eta) * param.RS;

% collateral constraints
param.lambda = 3;

% scale parameter technology
param.Aprod = 1;
param.Bprod = 1.132;
param.BU    = 1.3;
param.BP    = 2.3;

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