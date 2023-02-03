function hjb = HJB(V, G, param)

num0 = 1e-8; % numerical 0 for upwind scheme

[VaF, VaB] = deal(zeros(G.J, param.discrete_types));
for j = 1:param.discrete_types
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

cF = param.u1inv(VaF);
cB = param.u1inv(VaB);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

% IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
% IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
% I0 = ~IF & ~IB;
IF = (sF > 0) .* (G.grid(:, 1) < 1);  
IB = (sB < 0) .* (G.grid(:, 1) > 0) .* (IF == 0);
I0 = (1 - IF - IB);

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c) - param.v(G.N);

hjb.c = c; hjb.s = s; hjb.u = u;

end

