function [parents, parent_levels] = add_parents(grid, levels)

    d = size(grid, 2);
    maxl = max(max(levels)) + 1;
    [parent_mat, full_grid_vec, l_vec] = find_parents(maxl);
    
    parents_cell = cell(maxl * d, 1);
    parent_levels_cell = cell(maxl * d, 1);
    parents_old = grid;
    parent_levels_old = levels;
    for i = 1:(maxl * d)
        [parents, parent_levels] = add_same_dim_parents( ...
            parents_old, parent_levels_old, parent_mat, full_grid_vec, l_vec, maxl);
        % If new parent was already in parents_old, drop it
        % This command also removes duplicates in parents
        % At this point, parents_old has:
        % (1) been used in add_same_dim_parents, so don't need to again
        % (2) been added to parents_cell, so don't need to add
        [parents, IA] = setdiff(parents, parents_old, 'rows');
        parent_levels = parent_levels(IA, :);
        
        if isempty(parents)
            break;
        end
        
        % Delete duplicate points generated above
%         [parents, IA] = unique(parents, 'rows');
%         parent_levels = parent_levels(IA, :);
        parents_cell{i} = parents;
        parent_levels_cell{i} = parent_levels;
        
        parents_old = parents;
        parent_levels_old = parent_levels;
    end
    parents = cell2mat(parents_cell);
    parent_levels = cell2mat(parent_levels_cell);
end

function [parents, parent_levels] = ...
    add_same_dim_parents(grid, levels, parent_mat, full_grid_vec, l_vec, maxl)
    d = size(grid, 2);
    nmax = 2^maxl;
    
    parents_cell = cell(d, 1);
    parent_levels_cell = cell(d, 1);
    for k = 1:d
        idx = grid(:, k) * nmax + 1;
        % Find parents for each of the idx we just found
        [grid_idx, parent_idx] = find(parent_mat(idx, :));
        temp = grid(grid_idx, :);
        temp(:, k) = full_grid_vec(parent_idx);
        parents_cell{k} = temp;
        
        temp = levels(grid_idx, :);
        temp(:, k) = l_vec(parent_idx);
        parent_levels_cell{k} = temp;
    end
    
    parents = cell2mat(parents_cell);
    parent_levels = cell2mat(parent_levels_cell);
end

function [parent_mat, full_grid_vec, l_vec] = find_parents(maxl)
    nmax = 2^maxl;
    % Constructs a mapping of grid node (in 1D) to its level
    full_grid_cell = cellfun(@(l) (...
        0:2^-l:1)', num2cell(0:maxl)', 'UniformOutput', 0);
    % Manually code first two levels
    full_grid_cell{1} = 1/2; full_grid_cell{2} = [0; 1/2; 1];
    l_cell = cell(maxl + 1, 1);
    l_cell{1} = 0;
    for i = maxl+1:-1:2
        full_grid_cell{i} = setdiff(full_grid_cell{i}, full_grid_cell{i-1});
        l_cell{i} = repmat(i-1, numel(full_grid_cell{i}), 1);
    end
    full_grid_vec = cell2mat(full_grid_cell);
    l_vec = cell2mat(l_cell);
    [full_grid_vec, I] = sort(full_grid_vec);
    l_vec = l_vec(I);
    
    % Find all ancestors (parents) for each node
    % (i, j) in parent_mat denotes that node j is parent of node i
    i_cell = cell(maxl, 1);
    j_cell = cell(maxl, 1);
    i_cell{1} = [0; 1]; % Computes parent for level 1 manually
    j_cell{1} = [0.5; 0.5];
    for l = 2:maxl
        h = 2^-l;
        temp = full_grid_cell{l+1};
        i_cell{l} = repmat(temp, 2, 1);
        j_cell{l} = [temp + h; temp - h];
    end
    parent_mat = sparse(cell2mat(i_cell) * nmax + 1, ...
        cell2mat(j_cell) * nmax + 1, 1, nmax + 1, nmax + 1);
    for l = 1:maxl
        parent_mat = parent_mat * (parent_mat + speye(nmax + 1));
    end
    parent_mat = min(parent_mat, 1);
    
end