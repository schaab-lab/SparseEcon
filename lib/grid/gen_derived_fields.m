function G = gen_derived_fields(G)

    G.h = 2.^-G.lvl;
    G.J = size(G.grid, 1);
    [~, G.H_comp] = gen_H_mat(G.grid, G.lvl);

    % Economic units, defined by min and max inputs
    G.grid_to_value = @(x) x .* G.range + G.min;
    G.value = G.grid_to_value(G.grid);

    % For sufficiently small grids, use full matrices instead
    G.sparse = (G.J > 100);

    if isfield(G, 'names') && isfield(G, 'named_dims')
        for i = 1:numel(G.names)
            G.(G.names{i}) = G.value(:, G.named_dims{i});
            G.(['d', G.names{i}]) = ...
                G.range(G.named_dims{i}) .* min(G.h(:, G.named_dims{i}));
        end
    end

    G = gen_bound_grid(G);
    G = gen_FD_interior(G);

end
