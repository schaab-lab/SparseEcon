function G = gen_FD_interior(G)

    % Computing finite difference operators
    grid_to_bound = G.bound_grid_cell.grid_to_bound;

    for k = 1:G.d
        bound_grid = G.bound_grid_cell.grid{k};
        bound_lvl = G.bound_grid_cell.lvl{k};
        BH_grid_to_bound_comp = G.bound_grid_cell.BH_comp(G.bound_grid_cell.ids{k}, :);
        n_bound = size(bound_grid, 1);

        [left_neighbor_bound, right_neighbor_bound] = find_neighbors(bound_grid, 1:n_bound, k);
        interior_idx = find(bound_grid(:, k) > 0 & bound_grid(:, k) < 1);

        left_dist = bound_grid(1:n_bound, k) - bound_grid(left_neighbor_bound, k);
        right_dist = bound_grid(right_neighbor_bound, k) - bound_grid(1:n_bound, k);

        [~, bound_Hk] = gen_H_mat2(bound_grid, bound_lvl, k);
        bound_Ek = inv(bound_Hk);

        DF = sparse([interior_idx; interior_idx], ...
            [interior_idx; right_neighbor_bound(interior_idx)], ...
            [-1./right_dist(interior_idx), 1./right_dist(interior_idx)], ...
            n_bound, n_bound, 2*n_bound);
        D1F_interior = bound_Ek * DF * bound_Hk * BH_grid_to_bound_comp / G.range(k);
        G.DS_interior.D1F{k} = D1F_interior(grid_to_bound{k}, :);

        DB = sparse([interior_idx; interior_idx], ...
            [interior_idx; left_neighbor_bound(interior_idx)], ...
            [1./left_dist(interior_idx), -1./left_dist(interior_idx)], ...
            n_bound, n_bound, 2*n_bound);
        D1B_interior = bound_Ek * DB * bound_Hk * BH_grid_to_bound_comp / G.range(k);
        G.DS_interior.D1B{k} = D1B_interior(grid_to_bound{k}, :);

        [a, b, c] = stencil_central1(left_dist(interior_idx), right_dist(interior_idx));
        DC = sparse([interior_idx; interior_idx; interior_idx], ...
            [left_neighbor_bound(interior_idx); interior_idx; right_neighbor_bound(interior_idx)], ...
            [a; b; c], ...
            n_bound, n_bound, 3*n_bound);
        D1C_interior = bound_Ek * DC * bound_Hk * BH_grid_to_bound_comp / G.range(k);
        G.DS_interior.D1C{k} = D1C_interior(grid_to_bound{k}, :);

        if ismember(k, G.dxx_dims)
            [a, b, c] = stencil_central2(left_dist(interior_idx), right_dist(interior_idx));
            D2 = sparse([interior_idx; interior_idx; interior_idx], ...
                [left_neighbor_bound(interior_idx); interior_idx; right_neighbor_bound(interior_idx)], ...
                [a; b; c], ...
                n_bound, n_bound, 3*n_bound);
            D2_interior = bound_Ek * D2 * bound_Hk * BH_grid_to_bound_comp / G.range(k)^2;
            G.DS_interior.D2{k} = D2_interior(grid_to_bound{k}, :);
        end
    end

end