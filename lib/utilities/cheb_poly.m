function T = cheb_poly(n, x, a, b) 

    % x should be a 1x m vector at which the Chebyshev polynomials are to be evaluated

    x = 2*(x - a)/(b-a) -1;

    T = zeros(n, length(x));
    T(1, :) = 1; %Cheb Poly of degree1
    if n > 1
        T(2, :) = x; %Cheb Poly of degree2
    end

    if n > 2
        tm2 = T(1, :);
        tm1 = T(2, :);
        for i = 3:n
            T(i, :) = 2*x.*tm1 - tm2;
            tm2 = tm1;
            tm1 = T(i, :);
        end
    end

    %disp('ChebPoly Matrix is Poly# X Node# --- fnct approx should be: phi * T, for phi 1 X Poly#')
end