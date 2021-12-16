function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 10;
param.d = 1; param.d_idio = 1; param.d_agg = 0;

param.l_dense = 10; % vector of "surplus" for dense grid

param.nmin =  -1;
param.nmax = 500;

param.min = [param.nmin];
param.max = [param.nmax];

% Grid adaptation:
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 20;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_told\n'); end


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 500;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho = 0.05;
param.gamma = 2;
param.alpha = 0.33;
param.delta = 0.06;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

% Asset prices:
param.muR = 0.051;
param.r   = 0.030;

param.sigR  = sqrt(0.009722222222222);
param.sigR2 = 0.009722222222222;

% Earnings parameters:
param.zz  = [0.2, 1.5];
param.la1 = 1/3;
param.la2 = 1/3;
param.L   = param.la2/(param.la1+param.la2) * param.zz(1) + param.la1/(param.la1+param.la2) * param.zz(2);

param.discrete_types = numel(param.zz);


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