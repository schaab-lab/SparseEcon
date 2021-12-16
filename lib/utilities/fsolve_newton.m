function [x, J, obj, iter, step_size] = ...
    fsolve_newton(f, x0, obj0, y0, J0, max_counter_jacobian, display)

% Solver settings
tol_obj = 1e-9;
tol_dx = 1e-9;
maxit = 100;
step_size0 = 0.25;

if nargin <= 6, display = 2; end

counter_jacobian = 0;

if display >= 1
    fprintf('Initializing Newton solver. Obj tolerance: %g, dx tolerance: %g, maxit: %i.\n', ...
        tol_obj, tol_dx, maxit)
end
if ~exist('obj0', 'var')
    obj0 = f(x0, y0);
end
if any(isnan(obj0(:))), error('NaNs encountered at starting point.'), end
norm_obj = norm(obj0) / numel(obj0);

f_count = 0;
if exist('J0', 'var') && all(size(J0) == [numel(x0), numel(x0)])
    J = J0;
else
    time_jac = tic;
    J = compute_jacobian(f, x0, obj0, y0);
    time_jac = toc(time_jac);
    if display >= 1, fprintf('Time to compute initial Jacobian matrix: %.2f seconds.\n', time_jac); end
    f_count = f_count + numel(x0);
    counter_jacobian = counter_jacobian + 1;
end
if any(isnan(J(:))), error('Complex values in initial Jacobian construction.'), end

x = x0;
obj = obj0;
y = y0;
iter = 0;
step_size = step_size0;
dx_norm = tol_dx * 10;
last_iter_compute_jacobian = 0;
while (norm_obj > tol_obj) && (iter < maxit) && (dx_norm > tol_dx)
    norm_obj_old = norm_obj;
    x_old = x;
    y_old = y;
    obj_old = obj;
    dx = (-J\obj);
    dx_norm = norm(step_size*dx);
    x = x_old + step_size*dx;
    if numel(y0) == 4, [obj, y{1}, y{2}, y{3}, y{4}] = f(x, y); end
    if numel(y0) == 3, [obj, y{1}, y{2}, y{3}] = f(x, y); end
    if numel(y0) == 2, [obj, y{1}, y{2}] = f(x, y); end
    if numel(y0) == 1, [obj, y] = f(x, y); end
    norm_obj = norm(obj) / numel(obj);
    f_count = f_count + 1;
    
    % Requires at least tolObj of progress, to prevent an additional
    % iteration where objective function becomes nonresponsive
    % (Now, it will just keep decreasing step size until dxNorm -> 0)
    if any(isnan(obj)) || norm_obj > (norm_obj_old - tol_obj)
        step_size = 0.5*step_size;
        norm_obj = norm_obj_old;
        obj = obj_old;
        x = x_old;
        y = y_old;
        
        % If resulting dxNorm becomes too small, try recomputing J (if haven't already this iter)
        if ((dx_norm <= 100 * tol_dx) && ...
            (iter ~= last_iter_compute_jacobian || iter == 0) && ...
            (norm_obj > 10 * tol_obj))
            if display == 2, fprintf('Solver failed to make progress, recomputing Jacobian ...\n'); end
            J = compute_jacobian(f, x, obj, y);
            if any(isnan(J(:))), fprintf('Complex values in Jacobian reconstruction.\n'); break, end
            f_count = f_count + numel(x0);
            counter_jacobian = counter_jacobian + 1;
            step_size = step_size0;
            dx_norm = tol_dx * 10; % Force at least one more iteration
            last_iter_compute_jacobian = iter;
        end
    else
        iter = iter+1;
        if display == 2, fprintf('Iteration %.3i, fCount %i:  Objective = %9.5g,  stepSize = %9.5g,  norm(dx) = %9.5g.\n', ...
            iter, f_count, norm_obj, step_size, dx_norm); end
        
        delta = obj-obj_old;
        
        J = J + (delta/step_size - J*dx)*dx' / (dx'*dx);
        step_size = min(2*step_size, 2); %1.2*stepSize;
    end
    
    if nargin >= 6
        if counter_jacobian >= max_counter_jacobian
            fprintf('Terminating solver: Jacobian recomputed %.i times \n', max_counter_jacobian);
            break;
        end
    end
end
if (norm_obj <= tol_obj)
    if display >= 1, fprintf('Solver converged.\n\n'); end
elseif (iter >= maxit)
    if display >= 1, fprintf('Solver stopped after reaching maximum iteration allowance (%i).\n\n', maxit); end
elseif (dx_norm <= tol_dx)
    if display >= 1, fprintf(['Solver stopped after failing to make progress. Local minimum possible.\n', ...
        'Objective = %9.5g,  stepSize = %9.5g,  norm(dx) = %9.5g.\n\n'], ...
        norm_obj, step_size, dx_norm); end
else
    fprintf('Solver stopped for unknown reason.\n\n') 
end

end

function J = compute_jacobian(f, x0, obj0, y0)
    phi = 100;
    xi  = 0;
    
    h = phi * sqrt(eps) * max(abs(x0(:)), 1);
    J = zeros(numel(obj0), numel(x0));
    parfor k = 1:numel(x0)
        x1 = x0;
        x1(k) = x1(k) + h(k);
        obj1 = f(x1, y0);

        delta = obj1(:) - obj0(:);
        delta( abs(delta) < xi ) = 0;
        J(:, k) = delta / h(k);
    end
end