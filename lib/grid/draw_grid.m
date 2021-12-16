function draw_grid(grid, x, y, operator)

    colors = {'b', 'r', 'g', 'm', 'c', 'y'};
    offset = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01];

    assert(size(grid, 2) == 2, 'Error: Can only draw two-dimensional grids.')

    figure('Units', 'normalized', 'Position', [0, 1, 0.5, 0.5])
    hold on
    xlim([min(grid(:, 1))-0.05, max(grid(:, 1))+0.05]);
    ylim([min(grid(:, 2))-0.05, max(grid(:, 2))+0.05]);

    for j = 1:numel(x)
        i = find(grid(:, 1) == x(j) & grid(:, 2) == y(j));
        ids = find(operator(i, :));

        scatter(grid(ids, 1) + offset(j), grid(ids, 2), 100, 'filled', ...
            'MarkerFaceColor', colors{j})
        scatter(x(j), y(j), 300, 'o', 'MarkerEdgeColor', colors{j})
        scatter(x(j), y(j), 300, 'x', 'MarkerEdgeColor', colors{j})
        text(grid(ids, 1) + 2*offset(j), grid(ids, 2) + 2*offset(j), ...
            string(full(operator(i, ids))), 'FontSize', 14)
    end

    scatter(grid(:, 1), grid(:, 2), 'filled', 'MarkerFaceColor', 'k')
    xlabel('x'); ylabel('y')
    hold off

end