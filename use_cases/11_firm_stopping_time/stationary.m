function [diff, G, GDense, ss] = stationary(x, G, GDense, param)

%% AGGREGATES
r = x(1); if r > param.rho || r < -0.1, diff = NaN(1); return; end

Y = param.L;
w = 1;


%% VFI
G.income = r * G.a + w .* param.zz;

% State-constrained boundary conditions:
left_bound  = param.u1(G.income(G.grid(:,1)==0, :));
right_bound = param.u1(G.income(G.grid(:,1)==1, :));
for j = 1:param.discreteTypes
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
    GDense = gen_FD(GDense, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G,'V0'), G.V0 = param.u(G.income) / param.rho; end

% Solve VFI:
[V, hjb] = VFI(G, [], param);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
for j = 1:param.discreteTypes, muDense{j}  = G.BHDense * hjb.mu{j}; end
for j = 1:param.discreteTypes, sigDense{j} = G.BHDense * hjb.sig{j}; end

g = KF(muDense, sigDense, GDense, param);


%% MARKET CLEARING
B = sum(sum( GDense.a .* g .* GDense.dx));
C = sum(sum( (G.BHDense * hjb.c) .* g .* GDense.dx));
S = sum(sum( (G.BHDense * hjb.s) .* g .* GDense.dx));

diff = B;

ss.V = V; ss.g = g; ss.c = hjb.c; ss.s = hjb.s;
ss.B = B; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; ss.excess_supply = Y - C;
end




