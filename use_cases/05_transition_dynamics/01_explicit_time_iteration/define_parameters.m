function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 5;
param.d = 1; param.d_idio = 1; param.d_agg = 0;

param.l_dense = 8; % vector of "surplus" for dense grid

param.kmin = 0;
param.kmax = 100;

param.min = [param.kmin];
param.max = [param.kmax];

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


%% TRANSITION DYNAMICS PARAMETERS
param.time_grid_adjustment = 1;
param.T = 100; 
param.N = 120;
param.t = linspace(0, param.T, param.N)';
%{
    Here is what works: T=100, N=200, "nodal".
    For "cheb", here is what works: T = 100; NN = 100 / 200; H >= 25; it seems "cheb" does not like sparser time grids
    In some cases I need to tune H a little for "cheb" (what usually works: 20 - 25)
%}

param.bfun_type = "cheb"; 
param.cheb_H = 25;

if param.time_grid_adjustment == 1
    if param.N / param.T >= 2
        adjustment = @(x) x; 
    elseif param.N / param.T >= 1
        adjustment = @(x) (exp(x/param.T)-1) * param.T / (exp(1)-1);
    elseif param.N / param.T > 0.8
        adjustment = @(x) x.^2 / param.T^1;
    else 
        adjustment = @(x) x.^3 / param.T^2;
    end
    param.t = adjustment(param.t);
end
param.dt = diff(param.t); param.dt(param.N) = param.dt(param.N-1);

param.H(1) = param.N; if param.bfun_type == "cheb", param.H(1) = param.cheb_H; end
param.H(2) = 2; % # of time series to guess 

if param.N <= 6*param.T, param.implicit_g = 1; else param.implicit_g = 0; end


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho = 0.02;
param.gamma = 2;
param.eta   = 2;
param.alpha = 0.33;
param.delta = 0.025;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

param.v     = @(x) x.^(1+param.eta) / (1+param.eta);
param.v1    = @(x) x.^param.eta;
param.v1inv = @(x) x.^(1/param.eta);

% Earnings parameters:
param.zz  = [0.8, 1.2];
param.la1 = 1/3;
param.la2 = 1/3;
param.L   = param.la1/(param.la1+param.la2) * param.zz(1) + param.la2/(param.la1+param.la2) * param.zz(2);

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