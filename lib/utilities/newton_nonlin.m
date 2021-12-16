function [x, it] = newton_nonlin(f, fprime, x0, a, crit)

    x = x0; dx = 1*size(x0); it = 1;
    while (max(max(dx)) > crit || max(max(abs(f(x, a)))) > crit)  
        xnew = x - (f(x, a) ./ fprime(x));  
        dx = abs(x-xnew);
        x = xnew;
        it = it+1;
    end

end