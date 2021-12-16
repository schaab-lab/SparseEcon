function BH_comp = get_projection_matrix(points, lvl_points, G)
% Constructs matrix for projecting onto designated points
% INPUTS:
% - points: (m x d) matrix of points for projection
% - lvl_points: (m x d) matrix of levels associated with points
%               (optional - enter empty [] if points are not associated with grid)
% - G: Grid struct
%
% OUTPUTS:
% - In BH_comp, each row is a point, each column is a basis node
%   uH = BH_comp(i, :)*f gives the projected value of function f at point i
%

    if isempty(lvl_points)
        lvl_points = 1000 * ones(size(points));
    end

    % Often times, most (or all) of points is contained within grid, 
    % so some rows in H_basis can be constructed via shortcut:
    % projection is simply identity mapping for overlapping points
    [overlap_bool, overlap_grid_idx] = ismember(points, G.grid, 'rows');
    overlap_bound_idx = find(overlap_bool);
    BH_comp_overlap = sparse(overlap_bound_idx, ...
                             overlap_grid_idx(overlap_bool), ...
                             ones(nnz(overlap_bool), 1), ...
                             size(points, 1), ...
                             G.J);
    BH_comp = BH_comp_overlap;

    if nnz(~overlap_bool) > 0
        BH_comp_nonoverlap = H_basis(points(~overlap_bool, :), ...
                                     lvl_points(~overlap_bool, :), ...
                                     G.grid, G.lvl) * G.H_comp;
        BH_comp(~overlap_bool, :) = BH_comp_nonoverlap;
    end

end