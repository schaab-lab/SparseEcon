function test_setup_grid
    % 1. Output the field names
    % 2. Output grid, lvl, and economic units
    
    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)

    % 1D grid example
    G = setup_grid(4, 0, 0, 10);
    test_grid(G);

    % 2D grid example
    G = setup_grid(4, [2, 0], [0, 10], [10, 100], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    test_grid(G);
    assert(isequal([G.x, G.y], G.value));
end

function test_grid(G)
    disp(G)
    disp([G.grid, G.lvl, G.value])
end
