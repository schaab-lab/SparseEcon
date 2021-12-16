function [diff, G, G_dense, ss] = stationary(x, G, G_dense, param)

K = x(1);
% tau = x(2);


%% AGGREGATES
Y  = K.^param.alpha .* param.L.^(1-param.alpha);
w  = (1-param.alpha) .* Y ./ param.L;
rk = param.alpha .* Y ./ K;

r = rk - param.delta;


%% VFI
G.income = r * G.a + w .* param.zz; % + tau;

% State-constrained boundary conditions
left_bound  = param.u1(G.income);
% right_bound = param.u1(G.income) .* (G.t<param.tmax) + ...
%               param.u1(10000) .* (G.t==param.tmax & G.a>0) + param.u1(w.*param.zz).* (G.t==param.tmax & G.a==0);
right_bound = param.u1(G.income);

right_bound_T = 1e-8 * param.u(1e-5 + G.a) .* [1,1];

for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{2}.left.type = '0';   BC{2}.right.type = 'D';
    BC{1}.left.f  = @(points) sparse_project(left_bound(:,j),  points, G);
    BC{1}.right.f = @(points) sparse_project(right_bound(:,j), points, G);
    BC{2}.right.f = @(points) sparse_project(right_bound_T(:,j), points, G);
    G = gen_FD_dirichlet(G, BC, num2str(j));
    BC{2}.right.type = 'D_birth';
    G_dense = gen_FD_dirichlet(G_dense, BC, num2str(j));
end

% Initialize guess V0
if ~isfield(G,'V0')
    G.V0 = param.u(G.income);
    % G.V0(G.t == param.tmax, :) = right_bound_T(G.t == param.tmax, :);
    G.V0(G.t == param.tmax, :) = 1e-8 * param.u(1e-5 + G.a(G.t == param.tmax)) .* [1,1];
end

[V, hjb] = VFI(G, [], param);

% MAIN REASON FOR NON-CONVERGENCE SO FAR: What happens when you arrive at T
% with a=0? Numerically, seems like you need at least 1e-5 of buffer.
% That's because the 'D' BC doesn't allow you to "earn income" in the last
% "period".

%{
figure; hold on; 
plot(G.a(G.t == param.tmax), G.V0(G.t == param.tmax, 2)); 
plot(G.a(G.t == param.tmax), 1e-8 * param.u(1e-8 + G.a(G.t == param.tmax))); hold off;

figure; scatter3(G.a, G.t, hjb.c(:,1));
figure; scatter3(G.a, G.t, hjb.s(:,1));
figure; scatter3(G.a, G.t, hjb.c(:,2));
figure; scatter3(G.a, G.t, hjb.s(:,2));

%}

%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discrete_types, s_dense{j} = G.BH_dense * hjb.s(:,j); end

g = KF(s_dense, G_dense, param);


%% MARKET CLEARING
KS = sum(sum( G_dense.a .* g .* G_dense.dx));
C  = sum(sum( (G.BH_dense * hjb.c) .* g .* G_dense.dx));
S  = sum(sum( (G.BH_dense * hjb.s) .* g .* G_dense.dx));

excess_supply  = Y - C - param.delta * K;
excess_savings = S;
excess_capital = KS - K;

diff = excess_capital;

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.K = K; ss.w = w; 
ss.excess_supply = excess_supply; ss.excess_savings = excess_savings;
ss.excess_capital = excess_capital;

%% TESTING
%{
max(max(abs(  hjb.s - (G.income - hjb.c)  )))
max(max(abs(  hjb.s - (r * G.a + w .* param.zz - hjb.c)  )))
max(max(abs(  hjb.s -r*G.a - w.*param.zz + hjb.c  )))

sum(sum(   (hjb.s -r*G.a - w.*param.zz + hjb.c) .* g * GDense.dx  ))
sum(sum(hjb.s.*g*GDense.dx)) - sum(sum(r*G.a.*g*GDense.dx)) - sum(sum(w.*param.zz.*g*GDense.dx)) ...
    + sum(sum(hjb.c.*g*GDense.dx))

sum(sum(hjb.s.*g*GDense.dx)) - r*K - w.*param.L + sum(sum(hjb.c.*g*GDense.dx))

S - r*K - w.*param.L + C
S - (rk - param.delta)*K - w.*param.L + C
S - rk*K - w.*param.L + C + param.delta*K
S - Y + C + param.delta*K

sum(sum(hjb.s(G.t==param.tmax,:) .* g(G.t==param.tmax,:) * GDense.dx))
%}

end




