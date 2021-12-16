function test_adapt_grid
    % 1. Adapts grid and plots error rate
    % 2. Adapts grid with 'pct' rule

    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)
    
    add_tol = 10.^(-1.76);
    f = @(x) 1 ./ (abs(0.5 - x(:, 1).^4 - x(:, 2).^4) + 0.1);

    G_dense = setup_grid(13, [0, 0], [0, 0], [1, 1], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    G = setup_grid(4, [0, 0], [0, 0], [1, 1], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    
    error_dense = [];
    J_dense = [];
    f_dense = f(G_dense.value);

    BH = H_basis(G_dense.grid, G_dense.lvl, G.grid, G.lvl);
    f_hat = BH * G.H_comp * f(G.value);
    error = mean((f_dense - f_hat).^2);
    J = G.J;
    
    fprintf('Adaptation with the "tol" rule:\n')
    stats.n_change = 1;
    while stats.n_change > 0
        [G, ~, ~, stats] = adapt_grid(G, f(G.value), [], ...
            'AddTol', add_tol, 'KeepTol', -Inf);

        BH = H_basis(G_dense.grid, G_dense.lvl, G.grid, G.lvl);
        f_hat = BH * G.H_comp * f(G.value);
        error = [error; mean((f_dense - f_hat).^2)];
        J = [J; G.J];
    end

    error = error(J < 1500);
    J = J(J < 1500);

    for i = 4:11
        G_dense_i = setup_grid(i, [0, 0], [0, 0], [1, 1], ...
            'NamedDims', {1, 2}, 'Names', {'x', 'y'});

        BH = H_basis(G_dense.grid, G_dense.lvl, G_dense_i.grid, G_dense_i.lvl);
        f_hat = BH * G_dense_i.H_comp * f(G_dense_i.value);
        error_dense = [error_dense; mean((f_dense - f_hat).^2)];
        J_dense = [J_dense; G_dense_i.J];
    end

    %% Figure
    figure
    scatter3(G.x, G.y, f(G.value), 6, 'filled')
    view(45, 45);
    xlabel('x'); ylabel('y'); zlabel('f')
    grid on
    exportgraphics(gcf, './output/adaptation_f.eps');

    figure('visible', 'off');
    loglog(J, error, '-o', J_dense, error_dense, '-o', 'LineWidth', 2, 'MarkerSize', 4)
    legend({'Adaptive sparse grid', 'Dense grid'})
    grid on
    xlabel('Number of grid points')
    ylabel('MSE')
    exportgraphics(gcf, './output/adaptation_error.eps');
    
    %% Adaptation with 'pct' rule
    G = setup_grid(4, [0, 0], [0, 0], [1, 1], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    
    fprintf('\nAdaptation with the "pct" rule:\n')
    stats.n_change = 1;
    while stats.n_change > 0
        [G, ~, ~, stats] = adapt_grid(G, f(G.value), [], ...
            'AddTol', 5, 'KeepTol', 10.^-2, 'AddRule', 'pct');
    end
    
    %% Adaptation with 'num' rule
    G = setup_grid(4, [0, 0], [0, 0], [1, 1], ...
        'NamedDims', {1, 2}, 'Names', {'x', 'y'});
    
    fprintf('\nAdaptation with the "num" rule:\n')
    stats.n_change = 1;
    while stats.n_change > 0
        [G, ~, ~, stats] = adapt_grid(G, f(G.value), [], ...
            'AddTol', 10, 'KeepTol', 10.^-2, 'AddRule', 'pct');
    end
    
    disp(' ')
end
