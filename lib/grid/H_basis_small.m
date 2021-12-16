function phi = H_basis_small(point, grid, lvl_grid, h)
%
% Constructs H basis function FOR ONE POINT
%
%     if ~exist('h', 'var')
%         h = 2.^-lvl_grid;
%     end

    phi = max(1 - abs(point - grid) ./ h, 0);
    phi(lvl_grid == 0) = 1;
    phi = prod(phi, 2)';
    
end