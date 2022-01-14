function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.l = 0; param.surplus = [5,2];
param.d = 2; param.d_idio = 2; param.d_agg = 0;

param.l_dense = [7,4]; % vector of "surplus" for dense grid

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


%% TRANSITION DYNAMICS PARAMETERS
param.time_grid_adjustment = 1;
param.T = 100; 
param.N = 120;

param.bfun_type = "nodal"; 
param.cheb_H = 25;

param.H(1) = param.N; if param.bfun_type == "cheb", param.H(1) = param.cheb_H; end
param.H(2) = 1; % # of time series to guess 


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

% TFP shock:
param.shock_mean = 0;
param.shock_level = 0.01;
param.theta_shock = log(2);


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
param.t = linspace(0, param.T, param.N)';
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

if param.N <= 6*param.T, param.implicit_g = 1; else param.implicit_g = 0; end
param.H(1) = param.N; if param.bfun_type == "cheb", param.H(1) = param.cheb_H; end


end