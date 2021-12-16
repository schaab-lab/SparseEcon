function output = H_basis(points, lvl_points, grid, lvl_grid, exclude, idx_only)
% Constructs H basis function
%
% INPUTS:
% - points: (m x d) matrix of points for projection
% - lvl_points: (m x d) matrix of levels associated with points
%               (optional - enter empty [] if points are not associated with grid)
% - grid: (n x d) matrix of grid nodes
% - lvl_grid: (n x d) matrix of levels associated with grid
% - exclude: (n x 1) vector of grid nodes to assume Hbasis == 0 (optional)
% - idx_only: bool to return indices instead of sparse matrix
%
% OUTPUTS:
% - In BH, each row is a point, each column is a basis node
%   uH = BH(i, :)*aH gives the H representation of the function at point i
%

    n_points = size(points, 1);
    n_grid = size(grid, 1);
    if isempty(lvl_points)
        lvl_points = 1000 * ones(size(points));
    end
    if nargin >= 5
        grid = grid(~exclude, :);
        lvl_grid = lvl_grid(~exclude, :);
        new_idx_to_old = find(~exclude);
    end
    
    [~, lvl_points_groups, lvl_points_groups_idx, points_ic_sort] = ...
        split_by_level(points, lvl_points);
    [grid_groups, lvl_grid_groups, lvl_grid_groups_idx, grid_ic_sort] = ...
        split_by_level(grid, lvl_grid);
    output = cell(size(lvl_grid_groups, 1), 1);
    parfor i = 1:size(lvl_grid_groups, 1)
        lvl_temp = lvl_grid_groups(i, :);
        d_temp = lvl_temp > 0;
        grid_idx = lvl_grid_groups_idx{i};
        grid_temp = grid_groups{i}(:, d_temp);
        h_temp = 2.^-lvl_temp(:, d_temp);
        
        lvl_points_groups_temp = all(lvl_points_groups >= lvl_temp, 2);
        points_idx = points_ic_sort(cell2mat(lvl_points_groups_idx(lvl_points_groups_temp)'));
        points_temp = points(points_idx, d_temp);
        
        phi = ones(numel(points_idx), numel(grid_idx), nnz(d_temp));
        for k = 1:nnz(d_temp)
            phi(:, :, k) = ...
                max(1 - abs(points_temp(:, k) - grid_temp(:, k)') ./ h_temp(k), 0);
        end
        phi = prod(phi, 3);
        [x_idx, y_idx, phi] = find(phi);
        output{i} = [points_idx(x_idx(:)), grid_idx(y_idx(:)), phi(:)];
    end

    %% Putting results together
    output = cell2mat(output);
    output(:, 2) = grid_ic_sort(output(:, 2));
    % Convert back to original indices if some nodes on grid were excluded
    if nargin >= 5
        output(:, 2) = new_idx_to_old(output(:, 2));
    end
    if nargin < 6 || ~idx_only
        BH = sparse(output(:, 1), output(:, 2), output(:, 3), n_points, n_grid);
        output = BH;
    end
    
end

% Splits grid (or points) by grouping nodes by level
function [grid_groups, lvl_groups, group_idx, ic_sort] = split_by_level(grid, lvl)
    [lvl_lorted, IA] = sortrows(lvl);
    [~, ic_sort] = ismember(IA, 1:size(grid, 1));
    grid_sorted = grid(IA, :);
    [lvl_groups, ~, IC] = unique(lvl_lorted, 'rows');
    [~, first_idx_group] = ismember(1:size(lvl_groups, 1), IC);
    first_idx_group = first_idx_group - 1;
    obs_per_group = diff([first_idx_group, size(grid, 1)]);
%     assert(sum(obs_per_group) == size(grid, 1));
    
    group_idx = cell(1, size(lvl_groups, 1));
    grid_groups = cell(1, size(lvl_groups, 1));
    for i = 1:size(lvl_groups, 1)
        group_idx{i} = (1:obs_per_group(i))' + first_idx_group(i);
        grid_groups{i} = grid_sorted(group_idx{i}, :);
    end
end