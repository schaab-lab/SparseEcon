function H = get_jacobians_market_clearing(x0, z0, f, param, query, display)


%% SETUP and INITIALIZATION
obj0 = f(x0, z0);

hx = param.phi_jacobian * sqrt(eps) * max(abs(x0(:)), 1);
hz = param.phi_jacobian * sqrt(eps) * max(abs(z0(:)), 1);


%% H_x: [K x K]
if any(strcmp(query, 'H_x'))
    if display > 0, fprintf('Computing H_x  ...'); run_time = tic; end
    J = zeros(numel(obj0), numel(x0));
    parfor k = 1:numel(x0)
        x1    = x0;
        x1(k) = x1(k) + hx(k);
        obj1  = f(x1, z0);

        delta  = obj1(:) - obj0(:);
        J(:,k) = delta / hx(k);
    end
    H.H_x = J;
    if display > 0, run_time = toc(run_time); fprintf(' finished in %.2f seconds\n', run_time); end
end


%% H_z: [K x N]
if any(strcmp(query, 'H_z'))
    if display > 0, fprintf('Computing H_z  ...'); run_time = tic; end
    J = zeros(numel(obj0), numel(z0));
    parfor k = 1:numel(z0)
        z1    = z0;
        z1(k) = z1(k) + hz(k);
        obj1  = f(x0, z1);

        delta  = obj1(:) - obj0(:);
        J(:,k) = delta / hz(k);
    end
    H.H_z = J;
    if display > 0, run_time = toc(run_time); fprintf(' finished in %.2f seconds\n', run_time); end
end


end

