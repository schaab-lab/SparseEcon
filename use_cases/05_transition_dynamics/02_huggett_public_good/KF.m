function g = KF(s, G, param)

Az = FD_operator(G, param.theta_z * (param.zmean - G.z), param.sig_z*ones(G.J,1), 2);
Aa = FD_operator(G, s, zeros(G.J,1), 1);
AT = (Aa + Az)';

b = zeros(G.J, 1);

i_fix = 1;
b(i_fix) = 0.1;
row = [zeros(1, i_fix-1), 1, zeros(1, G.J-i_fix)];
AT(i_fix, :) = row;

gg = AT \ b;
g_sum = gg' * ones(G.J, 1) * G.dx;
g = gg ./ g_sum;

mass = sum(g * G.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end


end