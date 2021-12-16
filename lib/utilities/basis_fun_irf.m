function [X, nodes] = basis_fun_irf(Y, PHI, H, N, bfun_type, t, direction)
%{
    Function inputs:
        - Y: matrix of functions on which to project
        - PHI: A matrix of coefficients which to use to build functions
        - N: Number of functions
        - bfun_type: Type of basis
        - direction: what function will compute 
    Outputs:
        - Based on "direction": either matrix of functions or matrix of
          coefficients.

        - Matrix of functions is length(t) X N
        - Matrix of coefficients is 1 X H*N
%}

    % IMPORTANT:
    t = reshape(t, [1, length(t)]);

    % Construct vector of: nodes (= row vector)
    switch bfun_type    
        case "cheb"
            nodes = cheb_nodes(H, t(1), t(end)); 

        case "nodal"
            nodes = t;
    end


    % Construct matrix of: basis functions
    switch bfun_type    
        case "cheb"
            T_node = cheb_poly(H, nodes, t(1), t(end));
            T_fine = cheb_poly(H,     t, t(1), t(end));
        case "nodal"
            T_node = basis_function(H, nodes, t(1), t(end), bfun_type, nodes, 1);
            T_fine = basis_function(H,     t, t(1), t(end), bfun_type,     t, 1);
    end

    % Retrieve coefficients via projection
    if direction == "get_coefficient"
        X = zeros(1, H*N);

        for i = 1:N
            X(1, H*(i-1)+1 : i*H) = interp1(t, Y(:, i), nodes) / T_node;
        end
    end


    % Compute functions using basis functions
    if direction == "get_function"

        for i = 1:N
            X(:, i) = (PHI(1, H*(i-1)+1 : i*H) * T_fine)';
        end

    end

end