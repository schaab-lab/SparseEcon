function g = KF(s, G, param)

Az = [-speye(G.J)*param.la1,  speye(G.J)*param.la1; ...
       speye(G.J)*param.la2, -speye(G.J)*param.la2];

Aa{1} = FD_operator(G, s(:,1), zeros(G.J,1), 1, '1');
Aa{2} = FD_operator(G, s(:,2), zeros(G.J,1), 1, '2');

AT = (blkdiag(Aa{1}, Aa{2}) + Az)';

b = zeros(param.discrete_types * G.J, 1);

i_fix = 1;
b(i_fix) = 0.1;
row = [zeros(1, i_fix-1), 1, zeros(1, param.discrete_types*G.J-i_fix)];
AT(i_fix, :) = row;

gg = AT \ b;
g_sum = gg' * ones(param.discrete_types*G.J, 1) * G.dx;
gg = gg ./ g_sum;

g = [gg(1:G.J), gg(G.J+1:2*G.J)];

mass = sum(g * G.dx);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end


end