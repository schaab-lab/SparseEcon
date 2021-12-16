function [a, b, c] = stencil_central1(dx_left, dx_right)

    % From https://drive.google.com/file/d/0B81VL20ggLWyVEdLalF6R0NFdWM/view
    a = -dx_right ./ (dx_left .* (dx_left + dx_right));
    b = (dx_right - dx_left) ./ (dx_left .* dx_right);
    c = dx_left ./ (dx_right .* (dx_left + dx_right));
    
end