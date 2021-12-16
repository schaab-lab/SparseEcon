function [a, b, c] = stencil_central2(dx_left, dx_right)

    % From Sundqvist and Veronis (1970)
    % and https://drive.google.com/file/d/0B81VL20ggLWyVEdLalF6R0NFdWM/view
    % Second derivative operator
    a = 2 ./ (dx_left .* (dx_left + dx_right));
    b = -2 ./ (dx_left .* dx_right);
    c = 2 ./ (dx_right .* (dx_left + dx_right));
    
end