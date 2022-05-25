function param = define_parameters()

% Grid parameters
param.l = 5;
param.d = 1; param.d_idio = 1; param.d_agg = 0;

param.xmin = 0.1;
param.xmax = 1.0;

param.min = [param.xmin];
param.max = [param.xmax];

% ASG parameters
param.add_rule = 'tol';
param.add_tol = 1e-5;
param.keep_tol = 1e-6; 
param.max_adapt_iter = 20;
if param.keep_tol >= param.add_tol, error('keep_tol should be smaller than add_tol\n'); end

% PDE parameters
param.Delta = 1000;
param.maxit = 100;
param.crit  = 1e-8;

param.Delta_KF = 1000;
param.maxit_KF = 100;
param.crit_KF  = 1e-8;

% Economic parameters
param.rho = 0.05;
param.gamma = 0.5;

param.theta_x = -0.01;
param.sig_x = 0.01;

param.outside_option = @(x) 10 + 0*x;

param.u     = @(x) x.^(1-param.gamma);% / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

end