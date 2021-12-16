function test_keep_dims
    % 1. Assert that outputted relationships between old and new hold
    % 2. Assert correct transfer of 'dxx_dims' and 'dxy_dims'
    % 3. Output the field names
    % 4. Output grid, lvl, and economic units
    
    fprintf('===== Log output for %s =====\n', mfilename)
    rng(05292011)

    % 3D grid example
    d = 3;
    G = setup_grid(4, zeros(1, d), zeros(1, d), ones(1, d), ...
        'NamedDims', {1, [2, 3]}, 'Names', {'x', 'y'}, ...
        'DxxDims', [1], 'DxyDims', [1, 3; 2, 3]);
    
    dims_kept = 2:3;
    [G_new, old_to_new, new_to_old] = keep_dims(G, dims_kept);
    
    assert(isempty(G_new.dxx_dims))
    assert(isequal(G_new.dxy_dims, [1, 2]))
    assert(isequal(G_new.names, {'y'}))
    assert(isequal(G_new.named_dims, {[1, 2]}))
    
    assert(isequal(G.grid(old_to_new, dims_kept), G_new.grid))
    assert(isequal(G.lvl(old_to_new, dims_kept), G_new.lvl))
    assert(isequal(G.grid(:, dims_kept), G_new.grid(new_to_old, :)))
    assert(isequal(G.lvl(:, dims_kept), G_new.lvl(new_to_old, :)))
    
    disp(G_new)
    disp([G_new.grid, G_new.lvl, G_new.value])
end
