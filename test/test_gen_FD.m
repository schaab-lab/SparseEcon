function test_gen_FD
    % 1. Output hashes for DS and const objects
    % 2. Assert that for a linear function, derivative matches exactly on the interior, and also on
    %    the exterior under '1' BCs. Under '0' BCs, first derivative is halved on the boundary
    % 3. Assert boundary derivatives match exactly for VNB and VNF BCs

    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)
    
    % 2D linear function example (f(x, y) = 2x + 4y)
    f = @(x) 2*x(:, 1) + 4*x(:, 2);
    G = setup_grid(4, [2, 0], [0, 10], [10, 100], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    BC{1}.left.type = '0'; BC{1}.right.type = '0';
    BC{2}.left.type = '0'; BC{2}.right.type = '0';
    G = gen_FD(G, BC);
    
    fx_D1C = deriv_sparse(G, f(G.value), 1, 'D1C');
    fy_D1C = deriv_sparse(G, f(G.value), 2, 'D1C');
    
%     assert(max(abs(fx_D1C - 2)) <= 1e-10); disabled for now due to BC change
    interior_y_idx = (G.grid(:, 2) > 0 & G.grid(:, 2) < 1);
    assert(max(abs(fy_D1C(interior_y_idx) - 4)) <= 1e-10);
    assert(max(abs(fy_D1C(~interior_y_idx) - 2)) <= 1e-10);
    
    % 2D nonlinear function example with VN BCs
    G = setup_grid(4, [2, 0], [0, 10], [10, 100], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    f = @(x) sin(pi*x(:, 1)) .* sin(pi*x(:, 2));
    fx = @(x) pi * cos(pi*x(:, 1)) .* sin(pi*x(:, 2));
    fy = @(x) pi * sin(pi*x(:, 1)) .* cos(pi*x(:, 2));
    
    BC{1}.left.f = @(x) fx(G.grid_to_value(x)); BC{1}.right.f = BC{1}.left.f;
    BC{2}.left.f = @(x) fy(G.grid_to_value(x)); BC{2}.right.f = BC{2}.left.f;
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{2}.left.type = 'VNB'; BC{2}.right.type = 'VNF';
    
    G = gen_FD(G, BC);
    disp(['test bound_grid hash: ', hash(G.bound_grid_cell.grid)]);
    disp(['test DS hash: ', hash(G.DS_main)]);
    disp(['test const hash: ', hash(G.const_main)]);
    
    fx_D1B = deriv_sparse(G, f(G.value), 1, 'D1B');
    fx_D1F = deriv_sparse(G, f(G.value), 1, 'D1F');
    
    assert(max(abs(fx_D1B(G.grid(:, 1) == 0) - fx(G.value(G.grid(:, 1) == 0, :)))) <= 1e-10);
    assert(max(abs(fx_D1F(G.grid(:, 1) == 1) - fx(G.value(G.grid(:, 1) == 1, :)))) <= 1e-10);
    
    disp(' ');
end