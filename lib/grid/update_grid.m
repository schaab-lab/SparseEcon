function G = update_grid(G, grid, lvl)
% Update grid after adaptation
%
% INPUTS:
% - G: Grid struct
% - grid: New grid points
% - lvl: Levels associated with new grid points
%

    assert(size(grid, 2) == G.d, ...
        'Updated grid must have same dimensionality as previous grid.')
    G.grid = grid;
    G.lvl = lvl;

    G = gen_derived_fields(G);

end
