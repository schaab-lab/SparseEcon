function test_gen_H
    % 1. Output example HComp matrix
    % 2. Assert hierarchical coefficient at center of grid (0.5, 0.5) is equal to nodal coefficients

    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)

    % 1D grid example
    G = setup_grid(4, 0, 0, 10);
    test_grid(G, 'test_H1');

    % 2D grid example
    G = setup_grid(4, [2, 0], [0, 10], [10, 100], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    test_grid(G, 'test_H2');

    disp(' ');
end

function test_grid(G, filename)
    figure('visible', 'off');
    spy(G.H_comp)
    set(gcf, 'PaperUnits', 'normalized')
    set(gcf, 'PaperPosition', [0 0 1 1])
    exportgraphics(gcf, ['./output/', filename, '.eps'])
    disp([filename, ' hash: ', hash(G.H_comp)]);
    
    f = rand(G.J, 1);
    f_H = G.H_comp * f;
    center_node = all(G.grid == 0.5, 2);
    assert(all(f_H(center_node) == f(center_node)));
end