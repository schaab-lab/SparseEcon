function [H, H_comp] = gen_H_mat2(grid, lvl, k_omit)
% Inputs
%   - grid:  (n X d) matrix of the grid, where each row is coordinates
%   - level: (n X d) matrix of level of each node in grid
%
% Outputs
%   - H: H basis transformation matrix for each dimension
%

    n = size(grid, 1);
    d = size(grid, 2);
    H = cell(d, 1);
    H_comp = speye(n);

    for k = setdiff(1:d, k_omit)
        dim_min = min(grid(:, k));
        dim_max = max(grid(:, k));
        left = grid(:, k) == dim_min;
        right = grid(:, k) == dim_max;

        % When computing H, have to consider each level l separately
        H{k} = speye(n);
        % Level 0 is already done by speye command above
    for l = 1:max(lvl(:, k))
        % All grid points with "higher" level than l
        subgrid_idx = lvl(:, k) <= l;
        % CAN PROBABLY speed this up by only finding neighbors for current lvl
        [left_neighbor, right_neighbor] = find_neighbors(grid, subgrid_idx, k);

        % Note: Neighbors will be wrong for boundaries without following, but
        % they are never used for H and D matrix generation, so can be ignored
        left_neighbor(left) = NaN;
        right_neighbor(right) = NaN;

        % Only updates H for nodes in level l, since neighbor is level specific
        % and hence the neighbor vectors are only correct for current level
        % CAN MAKE this part faster by collecting indices then making matrix
        if l == 1
            idx_left = find(grid(:, k) == 0);
            idx_right = find(grid(:, k) == 1);
            H{k} = H{k} + sparse([idx_left; idx_right], ...
                [right_neighbor(idx_left); left_neighbor(idx_right)], -1, n, n);
        else
            indices = find(lvl(:, k) == l);
            H{k} = H{k} + sparse([indices; indices], ...
                [right_neighbor(indices); left_neighbor(indices)], -1/2, n, n);
        end
    end
        H_comp = H_comp * H{k};
    end
    
end