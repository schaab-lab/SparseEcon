function [G, BH_adapt, blacklist, stats] = adapt_grid(G, f, blacklist, varargin)
%{
adapt_grid.m: Updates grid based on adaptive sparse grid methods

INPUTS:
- G:         Sparse grid to adapt
- f:         Nodal values for function; if 'f' has more than one column 
             (e.g., multiple functions), the maximum hierarchical coefficient
             for each node is used for adaptation purposes.
- blacklist: Matrix of blacklisted nodes (see variable inputs for details)

VARIABLE INPUTS:
- AddRule: Rule for nodes chosen for grid expansion
  * 'tol': (Default) Expand nodes with hierarchical coefficient greater than 'AddTol'
  * 'pct': Expand 'AddTol' percent of nodes with greatest hierarchical coefficients
  * 'num': Expand 'AddTol' nodes with greatest hierarchical coefficients
  * NOTE: 'pct' and 'num' include nodes that may have been expanded before,
          so adaptation may eventually end if a sufficiently high number of
          nodes with large coefficients have been found.
- AddTol: Argument for 'AddRule' input
- KeepTol: Tolerance for grid coarsening - nodes with hierarchical
           coefficients smaller than 'KeepTol' will be dropped
- BlkListOpt: Rule for blacklist
  * 'children': (Default) Do not add children who are on the blacklist
  * 'parent':   Do not add children whose immediate parent(s) are on the blacklist
- AddDims: Adapt only in a subset of dimensions

%}

    %% Parse input parameters
    p = inputParser;
    addParameter(p, 'BlkListOpt', 'children')
    addParameter(p, 'AddTol', 1e-4, @isnumeric)
    addParameter(p, 'AddRule', 'tol')
    addParameter(p, 'KeepTol', 8e-5, @isnumeric)
    addParameter(p, 'AddDims', 1:size(G.grid, 2))
    parse(p, varargin{:});
    
    %% Begin grid adaptation
    grid = G.grid;
    lvl = G.lvl;
    aH = G.H_comp * f;
    
    if isequal(p.Results.AddRule, 'tol')
        add_idx = find(max(abs(aH), [], 2) > p.Results.AddTol*(max(f(:))-min(f(:))));
    elseif isequal(p.Results.AddRule, 'pct')
        [aH_top, add_idx] = maxk(max(abs(aH), [], 2), ...
            ceil(size(grid, 1) * p.Results.AddTol));
        add_idx(aH_top <= p.Results.KeepTol*(max(f(:))-min(f(:)))) = [];
    elseif isequal(p.Results.AddRule, 'num')
        [aH_top, add_idx] = maxk(max(abs(aH), [], 2), p.Results.AddTol);
        add_idx(aH_top <= p.Results.KeepTol*(max(f(:))-min(f(:)))) = [];
    end
    keep = max(abs(aH), [], 2) > p.Results.KeepTol*(max(f(:))-min(f(:)));
    blacklist = unique([blacklist; grid(~keep, :)], 'rows');
    stats.coarsened = numel(keep) - nnz(keep);
    stats.keep_tol = p.Results.KeepTol;
    
    fprintf('Evaluating %i points to refine on ...\n', numel(add_idx))
    
    % Find children of points we want to add
    stats.n_old = size(grid, 1);
    [new_points, new_levels] = add_points(grid, lvl, add_idx, p.Results.AddDims);
    % Delete points we want to coarsen
    grid_new = grid(keep, :);
    lvl_new = lvl(keep, :);

    if ~isempty(new_points)
        % Gather children we want to add that are not already on the coarsened grid
        [new_points, IA] = setdiff(new_points, grid_new, 'rows');
        new_levels = new_levels(IA, :);
        % Also remove points that are on the blacklist
        [new_points, IA] = setdiff(new_points, blacklist, 'rows');
        new_levels = new_levels(IA, :);
    
        % If specified, also remove points whose parents are on the blacklist
        if isequal(p.Results.BlkListOpt, 'parent')
            newH = 2.^-new_levels;
            has_parent = new_levels > 0;
            n_parents = 2 * sum(has_parent, 2);
            parents_cell = cell(size(new_points, 1), 1);
            parents_idx_cell = cell(size(new_points, 1), 1);
            for i = 1:size(new_points, 1)
                parents = repmat(new_points(i, :), n_parents(i), 1);
                parents(:, has_parent(i, :)) = parents(:, has_parent(i, :)) + ...
                    [diag(newH(i, has_parent(i, :))); -diag(newH(i, has_parent(i, :)))];
                parents_cell{i} = parents;
                parents_idx_cell{i} = i * ones(n_parents(i), 1);
            end
            parents = cell2mat(parents_cell);
            parents_idx = cell2mat(parents_idx_cell);
            temp = parents_idx(ismember(parents, blacklist, 'rows'));
            new_point_drop = unique(temp);
            new_points(new_point_drop, :) = [];
            new_levels(new_point_drop, :) = [];
    %     elseif isequal(p.Results.blacklist, 'all')
    %         for i = 1:size(new_points, 1)
    %             new_point_parents = add_parents(new_points(i, :), new_levels(i, :));
    %             temp = ismember(new_point_parents, blacklist, 'rows');
    %             new_point_keep(i) = any(temp);
    %         end
        end
    end
    
    stats.n_new_points = size(new_points, 1);

    grid_new = [grid_new; new_points];
    lvl_new = [lvl_new; new_levels];

    % Add parents of all points on new grid
    [parents, parent_levels] = add_parents(grid_new, lvl_new);
    if ~isempty(parents)
        [parents, IA] = setdiff(parents, grid_new, 'rows');
        parent_levels = parent_levels(IA, :);
    end
    stats.n_new_parents = size(parents, 1);

    grid_new = [grid_new; parents];
    lvl_new = [lvl_new; parent_levels];
    stats.n_new = size(grid_new, 1);
    stats.n_change = stats.n_new - stats.n_old;

    [grid_new, IA] = sortrows(grid_new);
    lvl_new = lvl_new(IA, :);
    
    BH_adapt = H_basis(grid_new, lvl_new, G.grid, G.lvl) * G.H_comp;
    G = update_grid(G, grid_new, lvl_new);
    
    %% Report summary of results
    fprintf('Mean level of abs(aH): %g. Maximum level of grid: %i. \n', ...
        mean(abs(aH(:))), max(max(lvl_new)))
    fprintf(['Nodes coarsened: %i, children added: %i, parents added: %i, \n' ...
             'Net change: %i -> %i = << %i >>, blacklisted: %i. \n'], ...
        stats.coarsened, stats.n_new_points, stats.n_new_parents, ...
        stats.n_old, stats.n_new, stats.n_change, size(blacklist, 1))
    
end