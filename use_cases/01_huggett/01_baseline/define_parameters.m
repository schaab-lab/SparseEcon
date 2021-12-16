function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 2; param.surplus = [2, 0];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [7, 4]; % vector of "surplus" for dense grid

param.amin = -1;
param.amax = 20;
param.zmin = 0.8;
param.zmax = 1.2;

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
param.rho = 0.02;
param.gamma = 2;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

% Earnings parameters:
param.discrete_types = 1; %numel(param.zz);
param.L = 1;

param.zmean = 1;
param.theta_z = 0.25;
param.sig_z = 0.01;


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