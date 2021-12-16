function deriv = deriv_sparse(G, f, k, operator, name)

if nargin < 5
    name = 'main';
end

if strcmp(operator, 'D2')
    assert(ismember(k, G.dxx_dims), ...
        sprintf(['Error: D2 operator requested in dimension %i. ', ...
        'dxx operator in dimension %i must be enabled. ', ...
        'See "DxxDims" optional input in setup_grid().'], k, k))
end

if strcmp(operator, 'D22')
    k1 = max(k);
    k2 = min(k);
    assert(~isempty(G.dxy_dims) && any(ismember([k1, k2; k2, k1], G.dxy_dims, 'rows')), ...
        sprintf(['Error: D22 operator requested in dimensions %i and %i. ', ...
        'dxy operator in dimensions %i and %i must be enabled. ', ...
        'See "DxyDims" optional input in setup_grid().'], ...
        k1, k2, k1, k2))
    deriv = G.(['DS_', name]).(operator){k1, k2} * f + G.(['const_', name]).(operator){k1, k2};
else
    deriv = (G.DS_interior.(operator){k} + G.(['DS_', name]).(operator){k}) * f + ...
        G.(['const_', name]).(operator){k};
end

end

