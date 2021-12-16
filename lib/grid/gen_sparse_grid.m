function [grid, lvl] = gen_sparse_grid(d, n, surplus)

    if nargin < 3
        surplus = zeros(1, d);
    end

    if d == 1
        [grid, lvl] = gen_sparse_grid(2, n, [surplus, 0]);
        grid = grid(lvl(:, 2) == 0, 1);
        lvl = lvl(lvl(:, 2) == 0, 1);
        return
    end

    l = n * ones(1, d) + surplus;
    % Generate the one-dimensional grid points for each possible level l
    % NOTE: grid_1D{i} stores the points for level l = i-1
    grid_1D = cellfun(@(l) 0:2^-l:1, num2cell(0:max(l)), 'UniformOutput', 0);
    grid_1D{1} = 1/2; grid_1D{2} = [0 1/2 1]; % Manually code first two levels
    for i = max(l)+1:-1:2
        grid_1D{i} = setdiff(grid_1D{i}, grid_1D{i-1});
    end
    grid_1D_size = cellfun(@numel, grid_1D);

    % For "normal" (non-surplus) dimensions, there can be a maximum of n
    % dimensions with nonzero levels. Thus we:
    % 1. Find all such possible combinations of dimensions. 
    % 2. Combine them with the surplus dimensions.
    % 3. For each composite combination, find the possible grid levels that
    %    satisfy the criterion for sparse grids (accounting for surpluses).
    %    These possible grid levels are the same across all composite combinations.
    % 4. Distribute these grid levels across all dimension combinations to
    %    produce the final set of levels across all dimensions.
    n_surplus = nnz(surplus);
    surplus_dims = find(surplus);
    normal_dims = find(~surplus);
    d_pos_level = min(n, d - n_surplus); % Max possible num of normal dims w/ >0 level
    normal_dim_combs = nchoosek(normal_dims, d_pos_level);

    surplus_levels = arrayfun(@(x) 0:(n+x), surplus(surplus_dims), 'UniformOutput', 0);
    normal_levels = repmat({0:n}, 1, d_pos_level);
    level_combs = ndgrid2([surplus_levels, normal_levels]);
    surplus_level_combs = level_combs(:, 1:n_surplus);
    normal_level_combs = level_combs(:, n_surplus+1 : end);

    % If we want to keep a full grid in some dimensions, ignore them for sparse
    % Brumm & Scheidegger's levels are +1 compared to ours, so our n is smaller
    % by one, and our |l| is smaller by d: |l| <= n+d-1 in our scheme is |l| <= n
    keep = sum(max(surplus_level_combs - surplus(surplus_dims), 0), 2) + ...
        sum(normal_level_combs, 2) <= n;
    surplus_level_combs = surplus_level_combs(keep, :);
    normal_level_combs = normal_level_combs(keep, :);
    levels_cell = cell(size(normal_dim_combs, 1), 1);
    for i = 1:size(normal_dim_combs, 1)
        temp = zeros(nnz(keep), d);
        temp(:, surplus_dims) = surplus_level_combs;
        temp(:, normal_dim_combs(i, :)) = normal_level_combs;
        levels_cell{i} = temp;
    end
    lvls = cell2mat(levels_cell);
    lvls = unique(lvls, 'rows'); % Some levels are duplicates (e.g., all 0's)
    lvl_size = grid_1D_size(lvls + 1);
    lvl_num_el = prod(lvl_size, 2);

    % For each combination of levels across dimensions, generate the grid
    % according to those levels; combine all such grids
    one_vec = ones(1, d); % pre-generate vector of ones
    grid_cell = cell(size(lvls, 1), 1);
    lvl_cell = cell(size(lvls, 1), 1);
    for i = 1:size(lvls, 1)
        grid_out = zeros(lvl_num_el(i), d);
        siz = lvl_size(i, :);
        for j = 1:d
            x = grid_1D{lvls(i, j) + 1};
            s = one_vec;
            s(j) = numel(x);
            x = reshape(x, s);
            s = siz;
            s(j) = 1;
            x = repmat(x, s);
            grid_out(:, j) = x(:);
        end
        grid_cell{i} = grid_out;
        lvl_cell{i} = repmat(lvls(i, :), size(grid_cell{i}, 1), 1);
    end
    grid = cell2mat(grid_cell);
    [grid, crosswalk] = sortrows(grid, d:-1:1);
    lvl = cell2mat(lvl_cell);
    lvl = lvl(crosswalk, :);

end
