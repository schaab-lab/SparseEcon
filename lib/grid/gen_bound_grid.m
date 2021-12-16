function G = gen_bound_grid(G)

    bound_ids = [];
    for k = 1:G.d
        new_nodes = [G.grid; G.grid];
        new_nodes(1:G.J, k) = 0;
        new_nodes(G.J+1:end, k) = 1;
        new_lvls = [G.lvl; G.lvl];
        new_lvls(:, k) = 1;
        bound_grid = [G.grid; new_nodes];
        bound_lvl = [G.lvl; new_lvls];

        [bound_grid, IA, IC] = unique(bound_grid, 'rows');
        bound_lvl = bound_lvl(IA, :);
        grid_to_bound = IC(1:G.J);

        bound_grid_cell.grid_to_bound{k} = grid_to_bound;
        bound_grid_cell.grid{k} = bound_grid;
        bound_grid_cell.lvl{k} = bound_lvl;
        bound_ids = [bound_ids; k*ones(size(bound_grid, 1), 1)];
    end

    bound_grid_all = cell2mat(bound_grid_cell.grid');
    bound_lvl_all = cell2mat(bound_grid_cell.lvl');
    [bound_grid_all, IA, IC] = unique(bound_grid_all, 'rows');
    bound_lvl_all = bound_lvl_all(IA, :);

    bound_grid_cell.ids = cell(G.d, 1);
    for k = 1:G.d
        bound_grid_cell.ids{k} = IC(bound_ids == k);
    end

    [bound_grid_cell.bound_all_Hk, bound_grid_cell.bound_all_H_comp] = ...
        gen_H_mat(bound_grid_all, bound_lvl_all);

    bound_grid_cell.BH_comp = ...
        get_projection_matrix(bound_grid_all, bound_lvl_all, G);

    bound_grid_cell.bound_grid_all = bound_grid_all;
    bound_grid_cell.bound_lvl_all = bound_lvl_all;
    bound_grid_cell.gridJ = G.J;

    G.bound_grid_cell = bound_grid_cell;

end