function [mat, bound_const] = FD_operator(G, mu, sigma, dims, BC_name)
% Constructs finite difference operators given grid, drift, and diffusion
%
% INPUTS:
% - G: Grid struct
% - mu: (J x d') matrix of drift terms
% - sigma: (J x d') matrix of diffusion terms
% - dims: (d' x 1) vector, subset of dimensions referenced by mu and sigma
% - BC_name: (Optional) name of boundary conditions to use
%
% OUTPUTS:
% - mat: FD operator matrix
% - bound_const: Residual constant terms
%

if nargin < 5
    BC_name = 'main';
end

% Shortcut for specific common case
if G.d == 1 && ~G.sparse && nargin == 2 && nargout == 1
    mat = sparse((mu .* (mu >= 0)) .* G.(['DFull_', BC_name]).D1F{1} + ...
                 (mu .* (mu <  0)) .* G.(['DFull_', BC_name]).D1B{1});
    return
end

if G.sparse
    mat = sparse(G.J, G.J);
else
    mat = zeros(G.J, G.J);
end

bound_const = zeros(G.J, 1);

%%%%%%%%%%%%%%%%%%    DRIFT    %%%%%%%%%%%%%%%%%%%
for i = 1:numel(dims)
    k = dims(i);
    if any(mu(:, i) ~= 0)
        if G.sparse
            ijs = [vec_x_spijs(mu(:,i) .* (mu(:,i) >= 0), G.(['DSijs_', BC_name]).D1F{k});
                   vec_x_spijs(mu(:,i) .* (mu(:,i) <  0), G.(['DSijs_', BC_name]).D1B{k})];
            mat_drift = sparse(ijs(:, 1), ijs(:, 2), ijs(:, 3), G.J, G.J);
        else
            mat_drift = (mu(:,i) .* (mu(:,i) >= 0)) .* G.(['DFull_', BC_name]).D1F{k} + ...
                        (mu(:,i) .* (mu(:,i) <  0)) .* G.(['DFull_', BC_name]).D1B{k};
        end
        mat = mat + mat_drift;
        
        if nargout > 1
            const_drift = mu(:,i) .* (mu(:,i) >= 0) .* G.(['const_', BC_name]).D1F{k} + ...
                          mu(:,i) .* (mu(:,i) <  0) .* G.(['const_', BC_name]).D1B{k};
            bound_const = bound_const + const_drift;
        end
    end
end

%%%%%%%%%%%%%%%%    DIFFUSION    %%%%%%%%%%%%%%%%%%%
for i = 1:numel(dims)
    k = dims(i);
    if any(sigma(:, i) ~= 0)
        assert(ismember(k, G.dxx_dims), ...
            sprintf(['Error: Nonzero SIGMA detected in dimension %i. ', ...
            'dxx operator in dimension %i must be enabled. ', ...
            'See "DxxDims" optional input in setup_grid().'], k, k))
        if G.sparse
            ijs = vec_x_spijs(1/2 * sigma(:,i).^2, G.(['DSijs_', BC_name]).D2{k});
            mat_diffusion = sparse(ijs(:, 1), ijs(:, 2), ijs(:, 3), G.J, G.J);
        else
            mat_diffusion = 1/2 * sigma(:,i).^2 .* G.(['DFull_', BC_name]).D2{k};
        end
        mat = mat + mat_diffusion;
        
        if nargout > 1
            const_diffusion = 1/2 * sigma(:,i).^2 .* G.(['const_', BC_name]).D2{k};
            bound_const = bound_const + const_diffusion;
        end
    end
end

%%%%%%%%%%%%%%%%   CROSS TERMS   %%%%%%%%%%%%%%%%%%%
for i = 1:numel(dims)
    k = dims(i);
    if any(sigma(:, i) ~= 0)
        for j = 1:i-1
            l = dims(j);
            assert(~isempty(G.dxy_dims) && any(ismember([k, l; l, k], G.dxy_dims, 'rows')), ...
                sprintf(['Error: Nonzero SIGMA detected in dimensions %i and %i. ', ...
                'dxy operator in dimensions %i and %i must be enabled. ', ...
                'See "DxyDims" optional input in setup_grid().'], ...
                k, l, k, l))
            if G.sparse
                ijs = vec_x_spijs(sigma(:,i) .* sigma(:,j), G.(['DSijs_', BC_name]).D22{k, l});
                mat_cross = sparse(ijs(:, 1), ijs(:, 2), ijs(:, 3), G.J, G.J);
            else
                mat_cross = (sigma(:,i) .* sigma(:,j)) .* G.(['DFull_', BC_name]).D22{k, l};
            end
            mat = mat + mat_cross;

            if nargout > 1
                const_cross = sigma(:,i) .* sigma(:,j) .* G.(['const_', BC_name]).D22{k, l};
                bound_const = bound_const + const_cross;
            end
        end
    end
end

end

