function G = gen_FD(G, BC, name)
% Constructs FD operators for given boundary conditions
%
% INPUTS:
% - G: Grid struct
% - BC: Boundary condition inputs. This is a (d x 1) cell, where each
%       element is a struct with the following properties
%     - left.type: type of left boundary condition
%     - right.type: type of right boundary condition
%     - left.f: values associated with left boundary condition
%     - right.f: values associated with right boundary condition
% - name: (Optional) name associated with boundary condition, useful when
%         using multiple boundary conditions for the same grid
%

if nargin < 3
    name = 'main';
end

% If BCs have not changed for a dimension and G.DS appears correctly populated, skip computation of that dimension
compute_dims = find(~arrayfun(@(k) isfield(G, ['BC_', name]) && ...
    isequal(G.(['BC_', name]){k}, BC{k}) && ...
    G.J == size(G.(['DS_', name]).D1F, k), 1:G.d));
G.(['BC_', name]) = BC;

BC_left = cellfun(@(x) x.left.type, G.(['BC_', name]), 'UniformOutput', 0);
BC_right = cellfun(@(x) x.right.type, G.(['BC_', name]), 'UniformOutput', 0);
% Only require f for BCs that are not '0' or '1'
f_bound_left = cell(G.d, 1);
f_bound_left(~strcmp(BC_left, '0') & ~strcmp(BC_left, '1')) = ...
    cellfun(@(x) x.left.f, G.(['BC_', name])(~strcmp(BC_left, '0') & ~strcmp(BC_left, '1')), ...
    'UniformOutput', 0);
f_bound_right = cell(G.d, 1);
f_bound_right(~strcmp(BC_right, '0') & ~strcmp(BC_right, '1')) = ...
    cellfun(@(x) x.right.f, G.(['BC_', name])(~strcmp(BC_right, '0') & ~strcmp(BC_right, '1')), ...
    'UniformOutput', 0);

%% Computing finite difference operators
grid_to_bound = G.bound_grid_cell.grid_to_bound;
const_value = cell(G.d, 1);

for k = compute_dims
    %% Setting up
    bound_grid = G.bound_grid_cell.grid{k};
    bound_lvl = G.bound_grid_cell.lvl{k};
    BH_grid_to_bound_comp = G.bound_grid_cell.BH_comp(G.bound_grid_cell.ids{k}, :);
    n_bound = size(bound_grid, 1);
    
    [left_neighbor_bound, right_neighbor_bound] = find_neighbors(bound_grid, 1:n_bound, k);
    left_idx_bound = find(bound_grid(:, k) == 0); n_left_bound = numel(left_idx_bound);
    right_idx_bound = find(bound_grid(:, k) == 1); n_right_bound = numel(right_idx_bound);
    boundary_idx = [left_idx_bound; right_idx_bound]; n_boundary = numel(boundary_idx);
    interior_idx = setdiff((1:n_bound)', [left_idx_bound; right_idx_bound]);
    
    left_dist = bound_grid(1:n_bound, k) - bound_grid(left_neighbor_bound, k);
    right_dist = bound_grid(right_neighbor_bound, k) - bound_grid(1:n_bound, k);
    
    % This is faster than multiplying all the other Hk matrices again
    bound_Hk = G.bound_grid_cell.bound_all_Hk{k}(G.bound_grid_cell.ids{k}, G.bound_grid_cell.ids{k}) \ ...
        G.bound_grid_cell.bound_all_H_comp(G.bound_grid_cell.ids{k}, G.bound_grid_cell.ids{k});
    bound_Ek = inv(bound_Hk);

    % Distance for exterior ghost nodes (not actually created in memory)
    left_offset = min(right_dist(left_idx_bound));
    right_offset = min(left_dist(right_idx_bound));
    left_dist(left_idx_bound) = left_offset;
    right_dist(right_idx_bound) = right_offset;
    [a1, b1, c1] = stencil_central1(left_dist, right_dist);
    [a2, b2, c2] = stencil_central2(left_dist, right_dist);

    %% Boundary operator construction
    % Forward difference
    DF_left = sparse([left_idx_bound; left_idx_bound], ...
        [left_idx_bound; right_neighbor_bound(left_idx_bound)], ...
        [-1./right_dist(left_idx_bound), 1./right_dist(left_idx_bound)], ...
        n_bound, n_bound, 2*n_bound);
    if ismember(BC_right{k}, {'0', 'VNF'})
        DF_right = sparse(n_bound, n_bound);
    end

    % Sparse operator
    DFH = bound_Ek * (DF_left + DF_right) * bound_Hk * BH_grid_to_bound_comp / G.range(k);
    G.(['DS_', name]).D1F{k} = DFH(grid_to_bound{k}, :);

    % Backward difference
    DB_right = sparse([right_idx_bound; right_idx_bound], ...
        [right_idx_bound; left_neighbor_bound(right_idx_bound)], ...
        [1./left_dist(right_idx_bound), -1./left_dist(right_idx_bound)], ...
        n_bound, n_bound, 2*n_bound);
    if ismember(BC_left{k}, {'0', 'VNB'})
        DB_left = sparse(n_bound, n_bound);
    end

    % Sparse operator
    DBH = bound_Ek * (DB_left + DB_right) * bound_Hk * BH_grid_to_bound_comp / G.range(k);
    G.(['DS_', name]).D1B{k} = DBH(grid_to_bound{k}, :);

    if ismember(BC_left{k}, {'0', 'VNB'})
        DC_left = sparse( ...
            [left_idx_bound; left_idx_bound; left_idx_bound], ...
            [left_idx_bound; left_idx_bound; right_neighbor_bound(left_idx_bound)], ...
            [a1(left_idx_bound), b1(left_idx_bound), c1(left_idx_bound)], ...
            n_bound, n_bound, 3*n_bound);
    end
    if ismember(BC_right{k}, {'0', 'VNF'})
        DC_right = sparse( ...
            [right_idx_bound; right_idx_bound; right_idx_bound], ...
            [left_neighbor_bound(right_idx_bound); right_idx_bound; right_idx_bound], ...
            [a1(right_idx_bound), b1(right_idx_bound), c1(right_idx_bound)], ...
            n_bound, n_bound, 3*n_bound);
    end
    DCH = bound_Ek * (DC_left + DC_right) * bound_Hk * BH_grid_to_bound_comp / G.range(k);
    G.(['DS_', name]).D1C{k} = DCH(grid_to_bound{k}, :);
    
    if ismember(k, G.dxx_dims)
        if ismember(BC_left{k}, {'0', 'VNB'})
            D2_left = sparse( ...
                [left_idx_bound; left_idx_bound; left_idx_bound], ...
                [left_idx_bound; left_idx_bound; right_neighbor_bound(left_idx_bound)], ...
                [a2(left_idx_bound), b2(left_idx_bound), c2(left_idx_bound)], ...
                n_bound, n_bound, 3*n_bound);
        end
        if ismember(BC_right{k}, {'0', 'VNF'})
            D2_right = sparse( ...
                [right_idx_bound; right_idx_bound; right_idx_bound], ...
                [left_neighbor_bound(right_idx_bound); right_idx_bound; right_idx_bound], ...
                [a2(right_idx_bound), b2(right_idx_bound), c2(right_idx_bound)], ...
                n_bound, n_bound, 3*n_bound);
        end
        D2H = bound_Ek * (D2_left + D2_right) * bound_Hk * BH_grid_to_bound_comp / G.range(k)^2;
        G.(['DS_', name]).D2{k} = D2H(grid_to_bound{k}, :);
    end
    
    %% Const term construction
    G.(['const_', name]).D1F{k} = zeros(n_bound, 1);
    G.(['const_', name]).D1B{k} = zeros(n_bound, 1);
    G.(['const_', name]).D1C{k} = zeros(n_bound, 1);
    G.(['const_', name]).D2{k} = zeros(n_bound, 1);
    G.(['DCH_', name]){k} = sparse([], [], [], n_bound, n_boundary, 2*n_boundary);
    
    const_c_left = zeros(numel(left_idx_bound), 1);
    const_c_right = zeros(numel(left_idx_bound), 1);
    
    % Note: No const terms for reflecting boundaries
    if strcmp(BC_left{k}, 'VNB')
        const_c_left = f_bound_left{k}(bound_grid(left_idx_bound, :)) * -left_offset;
        
        % VNB: fx(0) = (f(0) - f(-1)) / h => f(-1) = f(0) - h*fx(0)
        G.(['const_', name]).D1B{k}(left_idx_bound) = const_c_left * -1/left_offset;
        G.(['const_', name]).D1C{k}(left_idx_bound) = const_c_left .* a1(left_idx_bound);
        G.(['DCH_', name]){k}(left_idx_bound, 1:n_left_bound) = ...
            spdiags(a1, 0, n_left_bound, n_left_bound);
        if ismember(k, G.dxx_dims)
            G.(['const_', name]).D2{k}(left_idx_bound) = const_c_left .* a2(left_idx_bound);
        end
    end
    if strcmp(BC_right{k}, 'VNF')
        const_c_right = f_bound_right{k}(bound_grid(right_idx_bound, :)) * right_offset;
        
        % VNF: fx(J) = (f(J+1) - f(J)) / h => f(J+1) = f(J) + h*fx(J)
        G.(['const_', name]).D1F{k}(right_idx_bound) = const_c_right * 1/right_offset;
        G.(['const_', name]).D1C{k}(right_idx_bound) = const_c_right .* c1(right_idx_bound);
        G.(['DCH_', name]){k}(right_idx_bound, end-n_right_bound+1:end) = ...
            spdiags(c1, 0, n_right_bound, n_right_bound);
        if ismember(k, G.dxx_dims)
            G.(['const_', name]).D2{k}(right_idx_bound) = const_c_right .* c2(right_idx_bound);
        end
    end
    const_value{k} = [const_c_left; const_c_right];
    
    G.(['const_', name]).D1B{k} = G.(['const_', name]).D1B{k}(grid_to_bound{k});
    G.(['const_', name]).D1F{k} = G.(['const_', name]).D1F{k}(grid_to_bound{k});
    G.(['const_', name]).D1C{k} = G.(['const_', name]).D1C{k}(grid_to_bound{k});
    if ismember(k, G.dxx_dims)
        G.(['const_', name]).D2{k} = G.(['const_', name]).D2{k}(grid_to_bound{k});
    end

    G = create_aux_fields(G, k, name);
end

for k = 1:G.d
    % Cross-terms
    for k2 = 1:k
        if ~any(ismember([k, k2], compute_dims)) || isempty(G.dxy_dims) ...
                || ~any(ismember([k, k2; k2, k], G.dxy_dims, 'rows'))
            continue;
        end
        
        G.(['DS_', name]).D22{k, k2} = (G.DS_interior.D1C{k} + G.(['DS_', name]).D1C{k}) ...
            * (G.DS_interior.D1C{k2} + G.(['DS_', name]).D1C{k2});
        
        D22FC = (G.DS_interior.D1F{k} + G.(['DS_', name]).D1F{k}) * ...
            (G.DS_interior.D1C{k2} + G.(['DS_', name]).D1C{k2});
        D22BC = (G.DS_interior.D1B{k} + G.(['DS_', name]).D1B{k}) * ...
            (G.DS_interior.D1C{k2} + G.(['DS_', name]).D1C{k2});
        G.(['DS_', name]).D22{k, k2}(G.grid(:, k) == 0, :) = D22FC(G.grid(:, k) == 0, :);
        G.(['DS_', name]).D22{k, k2}(G.grid(:, k) == 1, :) = D22BC(G.grid(:, k) == 1, :);
        
        % Store DS matrices in sparse-index format for easy access later
        [i, j, s] = find(G.(['DS_', name]).D22{k, k2});
        G.(['DSijs_', name]).D22{k, k2} = [i, j, s];
        
        % For sufficiently small grids, use full matrices instead
        if ~G.sparse
            G.(['DFull_', name]).D22{k, k2} = full(G.(['DS_', name]).D22{k, k2});
        end
        
        % Shortcut for case when both boundaries are reflecting
        if strcmp(BC_left{k2}, '0') && strcmp(BC_right{k2}, '0')
            G.(['const_', name]).D22{k, k2} = zeros(G.J, 1);
            continue
        end
        % Each row of D1S.C{k} shows for that node, how much it loads on
        % other nodes on the grid. Each row of D1S.C{k2} * ghostC{k2} shows
        % for each node how much it needs add due to the ghost nodes.
        G.(['const_', name]).D22{k, k2} = (G.DS_interior.D1C{k} + G.(['DS_', name]).D1C{k}) * ...
            G.(['DCH_', name]){k2}(grid_to_bound{k2}, :) * const_value{k2};

        ghostD22FC = (G.DS_interior.D1F{k} + G.(['DS_', name]).D1F{k}) * ...
            G.(['DCH_', name]){k2}(grid_to_bound{k2}, :) * const_value{k2};
        ghostD22BC = (G.DS_interior.D1B{k} + G.(['DS_', name]).D1B{k}) * ...
            G.(['DCH_', name]){k2}(grid_to_bound{k2}, :) * const_value{k2};
        G.(['const_', name]).D22{k, k2}(G.grid(:, k) == 0) = ghostD22FC(G.grid(:, k) == 0);
        G.(['const_', name]).D22{k, k2}(G.grid(:, k) == 1) = ghostD22BC(G.grid(:, k) == 1);
    end
end

end

function G = create_aux_fields(G, k, name)
    if ismember(k, G.dxx_dims)
        fields = {'D1F', 'D1B', 'D1C', 'D2'};
    else
        fields = {'D1F', 'D1B', 'D1C'};
    end
    
    for field = fields
        field = field{:};
        % Store DS matrices in sparse-index format for easy access later
        [i_interior, j_interior, s_interior] = find(G.DS_interior.(field){k});
        [i_boundary, j_boundary, s_boundary] = find(G.(['DS_', name]).(field){k});
        G.(['DSijs_', name]).(field){k} = [[i_interior, j_interior, s_interior]; ...
            [i_boundary, j_boundary, s_boundary]];
        
        % For sufficiently small grids, use full matrices instead
        if ~G.sparse
            G.(['DFull_', name]).(field){k} = full(G.DS_interior.(field){k} + G.(['DS_', name]).(field){k});
        end
    end
end