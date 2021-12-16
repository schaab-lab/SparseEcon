function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 3; param.surplus = 0;
param.d = 1; param.d_idio = 1; param.d_agg = 0;

param.l_dense = 8; % vector of "surplus" for dense grid

param.kmin = 0;
param.kmax = 20;

param.min = param.kmin;
param.max = param.kmax;

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 15;
if param.keep_tol >= param.add_tol, error('keepTol should be smaller than addTold\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 100;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS
param.AH    = 0.6;
param.AL    = 0.4;
param.kappa = 2;
param.rho   = 0.05;
param.gamma = 2;
param.alpha = 1/3;
param.delta = 0.05;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);


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