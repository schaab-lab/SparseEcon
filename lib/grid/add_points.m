function [new_points, new_levels] = add_points(grid, levels, add_idx, dims)

    if ~exist('dims', 'var')
        dims = 1:size(grid, 2);
    end
    d = size(grid, 2);
    
    % Add children of points where adaptive expansion is needed
    new_points_cell = cell(size(add_idx, 1), 1);
    new_levels_cell = cell(size(add_idx, 1), 1);
    for i = 1:size(add_idx, 1)
        point = grid(add_idx(i), :);
        l = levels(add_idx(i), :);
        level_new = l + 1;
        h_new = 2.^-level_new;
        new_points = [point + diag(h_new); point - diag(h_new)];
        new_levels = repmat(l + eye(d), 2, 1);
        % Optional: Add only points in certain dimensions
        new_points = new_points([dims(:); d+dims(:)], :);
        new_levels = new_levels([dims(:); d+dims(:)], :);
        outside = min(new_points, [], 2) < 0 | max(new_points, [], 2) > 1;
        new_points_cell{i} = new_points(~outside, :);
        new_levels_cell{i} = new_levels(~outside, :);
    end
    new_points = cell2mat(new_points_cell);
    new_levels = cell2mat(new_levels_cell);
    [new_points, IA] = unique(new_points, 'rows');
    new_levels = new_levels(IA, :);
    
end