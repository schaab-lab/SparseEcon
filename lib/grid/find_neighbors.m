function [left_neighbor, right_neighbor] = find_neighbors(grid, subgrid_idx, k)

    n = size(grid, 1);
    dims = size(grid, 2);
    subgrid = grid(subgrid_idx, :);
    n_subgrid = size(subgrid, 1);
    % Sorts grid matrix by dim d, generating crosswalk and its inverse
    % Three sets of indices:
    %     - Original indices from full grid ('grid" object)
    %     - Subgrid indices
    %     - Sorted indices
    % We find neighbors of each element in sorted indices,
    % then convert to subgrid indices, then finally original indices
    original_to_subgrid = find(subgrid_idx);
    [~, subgrid_to_sorted] = sortrows(subgrid, [setdiff(1:dims, k), k]);
    [~, sorted_to_subgrid] = ismember(1:n_subgrid, subgrid_to_sorted);
    % After the sort, the "left" neighbor is just the immediately preceding 
    % element, the "right" neighbor is just the next element. But this is 
    % in sorted indices. Have to convert indices back to original indices.
    left_neighbor = zeros(n, 1);
    temp = [1, 1:(n_subgrid-1)]; % function from sorted -> sorted
    temp = subgrid_to_sorted(temp); % sorted -> subgrid
    temp = temp(sorted_to_subgrid); % subgrid -> subgrid
    temp = original_to_subgrid(temp); % subgrid -> original
    left_neighbor(subgrid_idx) = temp;
    
    right_neighbor = zeros(n, 1);
    temp = [2:n_subgrid, 1]; % function from sorted -> sorted
    temp = subgrid_to_sorted(temp); % sorted -> subgrid
    temp = temp(sorted_to_subgrid); % subgrid -> subgrid
    temp = original_to_subgrid(temp); % subgrid -> original
    right_neighbor(subgrid_idx) = temp;
    
end