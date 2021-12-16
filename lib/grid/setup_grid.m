function G = setup_grid(n, surplus, min, max, varargin)

%{
setup_grid.m: Initializes grid structure

INPUTS:
- n:       Level of sparse grid
- surplus: Dimension-specific level surplus
- min:     Dimension-specific minimums in economic units
- max:     Dimension-specific maximums in economic units

VARIABLE INPUTS:
- NamedDims: Cell of vectors of named dimensions
- Names:     Names for named dimensions
- DxxDims:   Vector of dimensions to compute dxx operators
- DxyDims:   (n x 2) matrix of dimensions to compute dxy operators

%}

    %% Parse input parameters
    p = inputParser;
    addParameter(p, 'NamedDims', {}, @iscell)
    addParameter(p, 'Names', {})
    addParameter(p, 'DxxDims', [])
    addParameter(p, 'DxyDims', [])
    parse(p, varargin{:});
    assert(numel(p.Results.NamedDims) == numel(p.Results.Names), ...
        'Named dimensions and names must have equal number of elements.')

    %% Set up grid
    reserved_names = {'min', 'max', 'range', 'grid', 'lvl', 'value', 'h', 'J', 'd', 'DS', 'DFull', ...
        'const', 'DFH', 'DBH', 'DCH', 'Names', 'NamedDims', 'H_comp', 'BH_dense'};
    assert(numel(min) == numel(max), 'Min and max must have equal number of elements')
    G.d = numel(min);

    G.min = min(:)'; G.max = max(:)'; G.range = G.max - G.min;
    if any(G.range == 0), warning('Grid min and max are equal for at least one dimension.'); end
    [G.grid, G.lvl] = gen_sparse_grid(G.d, n, surplus);

    if numel(p.Results.NamedDims) > 0
        for i = 1:numel(p.Results.NamedDims)
            assert(~any(strcmp(p.Results.Names{i}, reserved_names)), ...
                'Name "%s" is reserved.', p.Results.Names{i});
        end
        G.names = p.Results.Names;
        G.named_dims = p.Results.NamedDims;
    end

    G.dxx_dims = p.Results.DxxDims;
    G.dxy_dims = p.Results.DxyDims;

    G = gen_derived_fields(G);

end
