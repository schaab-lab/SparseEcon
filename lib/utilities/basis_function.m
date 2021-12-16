function [T, T_prime] = basis_function(n, nodes, region_min, region_max, type, X, Y)

if type == "nodal"
    % Here "X" is the vector of nodal basis points, "Y" is unused
    assert(issorted(X(:), 'strictascend'), ...
        'Nodal basis points must be sorted ascending and unique.');
    % Adds region min and max to the vector for griddedInterpolant
    % Note: These are not becoming basis points - it simply allows the
    % interpolant to work properly between the left- and right-most nodes
    % and the boundaries of the grid
    min_included = ismember(region_min, X);
    if ~min_included, X = [region_min; X(:)]; end
    max_included = ismember(region_max, X);
    if ~max_included, X = [X(:); region_max]; end
    T = zeros(n, numel(nodes));
    for ii = 1:n
        % Create "data" for each nodal basis point by making the associated "triangle"
        data = X(:) * 0;
        data(~min_included + ii) = 1;
        % The basis function is then just evaluating this triangle
        F = griddedInterpolant(X(:), data);
        T(ii, :) = F(nodes);
    end
    T_prime = 0;
end

if type == "nodal_g0"
    T = basis_function(n-1, nodes, region_min, region_max, 'nodal', X, Y);
    % Here "X" is the vector of nodal basis points. "Y" is a two-column
    % vector whose first column is the "x values" of g0, and whose second
    % column is the "y values" of g0.
    
    % NOTE: Since g0 is an additional basis function, the output will have
    % one more dimension than numel(X).
    F = griddedInterpolant(Y(:, 1), Y(:, 2)); % Interpolates supplied g0
    T = [T; F(nodes)];
    T_prime = 0;
end

if type == "cheb"
    T = cheb_poly(n, nodes, region_min, region_max); 
    T_prime = cheb_poly_diff(n, nodes, region_min, region_max, 1); 
end

if type == "cheb_g0"
    % Note "X" input is actually here a vector of x's and "Y" g0's
    assert(nargin == 7, 'Chebyshev basis with augmented g0 requires g0 inputs.')
    x = X;
    g0 = Y;
    g0_node = zeros(numel(nodes), size(g0, 2));
    assert(min(nodes) >= min(x) & max(nodes) <= max(x), ...
        'Supplied g0 does not contain nodes.')
    for ii = 1:numel(nodes)
        [closest_idx, dist] = knnsearch(x, nodes(ii), 'K', 2);
        g0_node(ii, :) = sum((sum(dist) - dist') .* g0(closest_idx, :)) ./ sum(dist);
    end
    T = (g0_node .* cheb_poly(n, nodes, region_min, region_max)')';
    
    %%% NEED TO INTERPOLATE HERE? 
%     [~, idx] = min(abs(x - repmat(nodes, size(x,1), 1)));
%     ID = find(ismember(GsOld(:,2:end), dd, 'rows'));
%     idx = arrayfun(@(y) find(x==y,1),nodes)
%     da = (max(x) - min(x)) / (numel(x)-1);
%     gaF_node = ( g0(ismember(x,nodes)+1) - g0(x==nodes) ) ./ da;
%     gaB_node = ( g0(x==nodes) - g0(x==nodes-1) ) ./ da;
%     ga_node = gaF_node .* (gaF_node>=0) + gaB_node .* (gaB_node<0);
    % ===== Begin Allen's new code =====
    ga_node = zeros(numel(nodes), size(g0, 2));
    for ii = 1:numel(nodes)
        g0_right = g0(x > nodes(ii), :);
        g0_left =  g0(x < nodes(ii), :);
        [right_idx, right_D] = knnsearch(x(x > nodes(ii)), nodes(ii), 'K', 1);
        [left_idx, left_D]   = knnsearch(x(x < nodes(ii)), nodes(ii), 'K', 1);
        % THIS TREATMENT OF BOUNDARIES IS VERY AD-HOC
        if isempty(left_D)
            ga_node(ii, :) = (g0_right(right_idx, :) - g0(x == nodes(ii), :)) / right_D;
        elseif isempty(right_D)
            ga_node(ii, :) = (g0(x == nodes(ii), :) - g0_left(left_idx, :)) / left_D;
        else
            ga_node(ii, :) = (g0_right(right_idx, :) - g0_left(left_idx, :)) / (right_D + left_D);
        end
    end
    % ===== End Allen's new code =====
    T_prime = (g0_node.*cheb_poly_diff(n, nodes, region_min, region_max, 1)' ...
              + ga_node.*cheb_poly(n, nodes, region_min, region_max)')';
end

if type == "gaussian"
    assert(nargin == 7, 'Gaussian basis specified, require X and Y inputs.')
    for i = 1:n
        T(i,:) = exp( -(nodes-X(i)).^2 / (2*Y(i)^2) );
        T_prime(i,:) = -(nodes-X(i))/Y(i)^2 .* T(i,:);
    end
end

end