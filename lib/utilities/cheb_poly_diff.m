function T = cheb_poly_diff(n, x, a, b, order) 

    % x should be a 1x m vector at which the Chebyshev polynomials are to be evaluated
    t = cheb_poly(n, x, a, b); %Needs to come before the x-transformation in the next line!

    %x = 2*(x - a)/(b-a) -1;

    T = zeros(n, length(x));
    T(1, :) = 0;
    if n > 0, T(2, :) = 2/(b-a); end

    if n > 1
        for i = 2:n-1
            T(i+1, :) = 4/(b-a)*t(i, :) + 2*(2*(x - a)/(b-a) -1).*T(i, :) - T(i-1, :);
        end
    end


    if order == 2, disp('Not yet written'); end
    % if order == 2
    %         TT = zeros(n,length(x));
    %         TT(1, :) = 0;
    %         if n>0, TT(2, :) = 0; end;
    % 
    %         if n > 1
    %             for i = 2:n-1
    %                 TT(i+1, :) = 4*T(i, :) + 2*x.*TT(i, :) - TT(i-1, :);
    %             end
    %         end
    %         T = TT;
    % end

end