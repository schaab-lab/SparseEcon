function g = KF(s, G, param)

Az = FD_operator(G, param.theta_z * (param.zmean - G.z), param.sig_z*ones(G.J,1), 2);
Aa = FD_operator(G, s, zeros(G.J,1), 1);
AT = (Aa + Az)';


% KF #1:
AT1 = AT;
b = zeros(G.J, 1);

i_fix = 1;
b(i_fix) = 0.1;
row = [zeros(1, i_fix-1), 1, zeros(1, G.J-i_fix)];
AT1(i_fix, :) = row;

gg = AT1 \ b;
g_sum = gg' * ones(G.J, 1) * G.dx;
g = gg ./ g_sum;


% KF #2:
% g1 = zeros(G.J, 1);
% g1(G.a == param.amin & G.z == param.zmin,:) = 1 / G.dx;
% 
% for n = 1:param.maxit_KF   
%     B = 1/param.Delta_KF .* speye(G.J) - AT;
%     b = g1 / param.DeltaKF;
% 
%     g_new = B\b;
% 
%     diff = max(abs( g1 - g_new ));
%     if diff < param.crit_KF, break; end
%     g1 = g_new;
% end
% if n == param.maxit_KF, fprintf('KF did not converge. Remaining Gap: %.2d\n',diff); end


% Some tests:
mass = sum(g .* G.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end
% if max(abs(g1 - g)) > 1e-5, fprintf('Distributions g1 and g2 do not align!\n'); end


end