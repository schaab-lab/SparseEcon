function [ijs] = vec_x_spijs(v, ijs)
    % Performs elementwise multiplication between a n-COLUMN vector and a
    % (n x n) sparse matrix (stored using the [i,j,s] format), and returns
    % a (n x 3) matrix of [i,j,s] which can be then fed into
    % sparse(i,j,s,n,n) to create the product matrix
    
    v_expand = v(ijs(:, 1));
    keep = v_expand ~= 0;
    ijs = ijs(keep, :);
    ijs(:, 3) = ijs(:, 3) .* v_expand(keep);
    
end