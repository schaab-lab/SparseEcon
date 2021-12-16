function test_parents
    % Interpolating a grid onto itself is an identity, but only if grid is properly set up with parents.
    % Test (a): Dropping random points from the grid destroys this identity IF parents aren't added back
    % Test (b): Adding parents back restores this identity

    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)
    d = 2;
    l = 8;
    drop_pct1 = 5;
    drop_pct2 = 50;

    % f = @(x, y) 1 ./ (abs(0.5 - x.^4 - y.^4) + 0.1);
    f = @(x, y) log(x + 1) .* log(y + 1);
    [grid, lvl] = gen_sparse_grid(d, l);
    v = f(grid(:, 1), grid(:, 2));

    scatter3(grid(:, 1), grid(:, 2), v);

    fprintf("Starting tests with grid of depth %i, of size %i.\n\n", l, size(grid, 1));

    %% Test case 1a: Drop random points from the grid
    fprintf("Test case 1a: Drop %i percent of random points from the grid.\n", drop_pct1);
    sample = randsample(size(grid, 1), ceil(size(grid, 1) * (1-drop_pct1/100)));
    grid_new = grid(sample, :);
    lvl_new = lvl(sample, :);
    v_new = f(grid_new(:, 1), grid_new(:, 2));

    [H, H_comp] = gen_H_mat(grid_new, lvl_new);
    aH = H_comp * v_new;
    BH = H_basis(grid_new, lvl_new, grid_new, lvl_new);
    v_interp = BH * aH;

    scatter3(grid_new(:, 1), grid_new(:, 2), v_interp);
    dist = abs(v_new - v_interp);
    errors = numel(find(dist >= 1e-10));
    avg_error = mean(dist(dist >= 1e-10));
    assert(errors > 0)
    fprintf(['Test case 1a completed. Grid size is %i. Interpolation errors: %i.' ...
        ' Average error: %4f.\n'], size(grid_new, 1), errors, avg_error);

    %% Test case 1b: Add parents for remaining points after test 1a
    fprintf("Test case 1b: Add parents for remaining points after test 1.\n");
    [parents, parent_levels] = add_parents(grid_new, lvl_new);
    grid_new = [grid_new; parents];
    lvl_new = [lvl_new; parent_levels];
    % Delete duplicate points generated above
    [grid_new, IA] = unique(grid_new, 'rows');
    lvl_new = lvl_new(IA, :);

    v_new = f(grid_new(:, 1), grid_new(:, 2));

    [~, H_comp] = gen_H_mat(grid_new, lvl_new);
    aH = H_comp * v_new;
    BH = H_basis(grid_new, lvl_new, grid_new, lvl_new);
    v_interp = BH * aH;

    assert(isequal(BH * H_comp, speye(size(BH * H_comp))));

    scatter3(grid_new(:, 1), grid_new(:, 2), v_interp);
    dist = abs(v_new - v_interp);
    errors = numel(find(dist >= 1e-10));
    avg_error = mean(dist(dist >= 1e-10));
    fprintf(['Test case 1b completed. Grid size is %i. Interpolation errors: %i.' ...
        ' Average error: %4f.\n\n'], size(grid_new, 1), errors, avg_error);

    %% Test case 2a: Drop random points from the grid
    fprintf("Test case 2a: Drop %i percent of random points from the grid.\n", drop_pct2);
    sample = randsample(size(grid, 1), ceil(size(grid, 1) * (1-drop_pct2/100)));
    grid_new = grid(sample, :);
    lvl_new = lvl(sample, :);
    v_new = f(grid_new(:, 1), grid_new(:, 2));

    [~, H_comp] = gen_H_mat(grid_new, lvl_new);
    aH = H_comp * v_new;
    BH = H_basis(grid_new, lvl_new, grid_new, lvl_new);
    v_interp = BH * aH;

    scatter3(grid_new(:, 1), grid_new(:, 2), v_interp);
    dist = abs(v_new - v_interp);
    errors = numel(find(dist >= 1e-10));
    avg_error = mean(dist(dist >= 1e-10));
    assert(errors > 0)
    fprintf(['Test case 2a completed. Grid size is %i. Interpolation errors: %i.' ...
        ' Average error: %4f.\n'], size(grid_new, 1), errors, avg_error);

    %% Test case 2b: Add parents for remaining points after test 2a
    fprintf("Test case 2b: Add parents for remaining points after test 1.\n");
    [parents, parent_levels] = add_parents(grid_new, lvl_new);
    grid_new = [grid_new; parents];
    lvl_new = [lvl_new; parent_levels];
    % Delete duplicate points generated above
    [grid_new, IA] = unique(grid_new, 'rows');
    lvl_new = lvl_new(IA, :);

    v_new = f(grid_new(:, 1), grid_new(:, 2));

    [~, H_comp] = gen_H_mat(grid_new, lvl_new);
    aH = H_comp * v_new;
    BH = H_basis(grid_new, lvl_new, grid_new, lvl_new);
    v_interp = BH * aH;

    assert(isequal(BH * H_comp, speye(size(BH * H_comp))));

    scatter3(grid_new(:, 1), grid_new(:, 2), v_interp);
    dist = abs(v_new - v_interp);
    errors = numel(find(dist >= 1e-10));
    avg_error = mean(dist(dist >= 1e-10));
    fprintf(['Test case 2b completed. Grid size is %i. Interpolation errors: %i.' ...
        ' Average error: %4f.\n'], size(grid_new, 1), errors, avg_error);
    
    disp(' ')
end
