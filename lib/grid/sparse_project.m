function x_project = sparse_project(x, points, G)

    I = size(points, 1);
    x_project = zeros(size(points,1), size(x,2));

    if I < 5000
        xH = G.H_comp * x;
        for i = 1:I        
            BH = H_basis_small(points(i,:), G.grid, G.lvl, G.h);
            x_project(i, :) = BH * xH;    
        end    
    else
        BH = get_projection_matrix(points, points+100, G);
        x_project = BH * x;
    end

end