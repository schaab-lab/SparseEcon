function test_gen_BH
    % 1. Output example BH matrix
    % 2. Assert projection of linear function has zero error
    % 3. Assert projecting a grid onto itself is an identity
    % 4. Assert projecting a grid onto a subset of itself is an identity on that subset
    % 5. Assert the lvl_points input can be made arbitrarily large with no change to output

    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)

    % 1D grid example
    G = setup_grid(4, 0, 0, 10);
    f = @(x) 2*x(:, 1);
    test_grid(G, f, 'test_BH1');

    % 2D grid example
    G = setup_grid(4, [2, 0], [0, 10], [10, 100], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    f = @(x) 2*x(:, 1) + 4*x(:, 2);
    test_grid(G, f, 'test_BH2');
    
    disp(' ');
end

function test_grid(G, f, filename)
    points_interp = rand(100, G.d);
    BH = H_basis(points_interp, 100 * ones(size(points_interp)), G.grid, G.lvl);
    figure('visible', 'off');
    spy(BH)
    set(gcf, 'PaperUnits', 'normalized')
    set(gcf, 'PaperPosition', [0 0 1 1])
    exportgraphics(gcf, ['./output/', filename, '.eps'])
    disp([filename, ' hash: ', hash(BH)]);
    
    assert(max(abs(BH * G.H_comp * f(G.grid) - f(points_interp))) < 1e-10);
    
    BH = H_basis(G.grid, G.lvl, G.grid, G.lvl);
    assert(isequal(BH * G.H_comp, speye(size(BH * G.H_comp))));
    
    % Using the first 10 elements here so that speye returns corresponding identity)
    BH = H_basis(G.grid(1:10, :), G.lvl(1:10, :), G.grid, G.lvl);
    assert(isequal(BH * G.H_comp, speye(size(BH * G.H_comp))));
    
    BH = H_basis(G.grid, 100 * G.lvl, G.grid, G.lvl);
    assert(isequal(BH * G.H_comp, speye(size(BH * G.H_comp))));
end