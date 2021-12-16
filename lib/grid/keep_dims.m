function [G_new, old_to_new, new_to_old] = keep_dims(G, dims)
% Keep a subset of dimensions on grid
% INPUTS:
% - G: Grid struct
% - dims: Dimensions to keep
%
% OUTPUTS:
% - G_new: New grid struct
% - old_to_new: ("IA") Mapping from old grid indices to new grid indices
% - new_to_old: ("IC") Mapping from new grid indices to old grid indices
% (See Matlab help file for unique() for mapping usage examples)
%

    G_new.d = numel(dims);
    G_new.min = G.min(dims);
    G_new.max = G.max(dims);
    G_new.range = G.range(dims);

    [G_new.grid, old_to_new, new_to_old] = unique(G.grid(:, dims), 'rows');
    G_new.lvl = G.lvl(old_to_new, dims);

    if isfield(G, 'names') && isfield(G, 'named_dims')
        % For named dimensions, only transfer if all dimensions in set have been kept
        G_new.named_dims = {};
        G_new.names = {};
        for i = 1:numel(G.names)
            [LIA, LOCB] = ismember(G.named_dims{i}, dims);
            if all(LIA)
                G_new.named_dims = [G_new.named_dims, {LOCB}];
                G_new.names = [G_new.names, G.names(i)];
            end
        end
    end

    % Mapping of old dimension numbering to new dimension numbering
    [~, new_dxx_dims] = ismember(G.dxx_dims, dims);
    new_dxx_dims(new_dxx_dims == 0) = [];

    [~, new_dxy_dims] = ismember(G.dxy_dims, dims);
    new_dxy_dims(any(new_dxy_dims == 0, 2), :) = [];

    G_new.dxx_dims = new_dxx_dims;
    G_new.dxy_dims = new_dxy_dims;

    G_new = gen_derived_fields(G_new);

end
