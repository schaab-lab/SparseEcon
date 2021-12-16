function BH_interp = agg_interp(grid, idx, grid_fine, idx_fine)

    [grid_agg, IA] = unique(grid(:, (size(grid_fine, 2) + 1) : end), 'rows');
    idx_agg = idx(IA, (size(grid_fine, 2) + 1) : end);
    grid_cross_cell = cell(size(grid_agg, 1), 1);
    idx_cross_cell = cell(size(grid_agg, 1), 1);
    J_fine = size(grid_fine, 1);
    for ii = 1:size(grid_agg, 1)
        grid_cross_cell{ii} = [grid_fine, repmat(grid_agg(ii, :), J_fine, 1)];
        idx_cross_cell{ii} = [idx_fine, repmat(idx_agg(ii, :), J_fine, 1)];
    end
    grid_cross = cell2mat(grid_cross_cell);
    idx_cross = cell2mat(idx_cross_cell);

    BH_interp = Hbasis(grid_cross, idx_cross, grid, idx);

end