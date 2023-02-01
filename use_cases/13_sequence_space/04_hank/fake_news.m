function [H, p] = fake_news(x0, z0, ss, G, G_dense, param)


%% SETUP
f = @(x, z) transition(x, z, ss, G, G_dense, param, 'fake_news');

sim0 = f(x0, z0);

hx = param.phi_jacobian * sqrt(eps) * max(abs(x0(:)), 1);
hz = param.phi_jacobian * sqrt(eps) * max(abs(z0(:)), 1);

N = param.N; K = numel(x0); J = param.discrete_types * G.J; dt = param.dt;

policies = {'s', 'c'};

assert(numel(z0) == N);
assert(numel(x0) == K);


%% JACOBIANS: POLICIES
zp = z0; zp(N) = zp(N) + hz(N); 
sim_zp = f(x0, zp);

sim_xp = cell(K/N, 1);
for h = 1:K/N
    xp = x0; xp(h*N) = xp(h*N) + hx(h*N);
    sim_xp{h} = f(xp, z0);
end

% POLICIES: c, s
for j = 1:numel(policies)
    for n = 1:N
        for k = 1:N
            if n > k
                p.([policies{j}, '_z']){n}(:, k) = zeros(J, 1);
                continue; 
            end
            p.([policies{j}, '_z']){n}(:, k) = ...
                (sim_zp.(policies{j}){N-(k-n)}(:) - sim0.(policies{j}){N-(k-n)}(:)) / hz(N);
        end
    end
end
for j = 1:numel(policies)
    for h = 1:K/N
        for n = 1:N
            for k = 1:N
                if n > k
                    p.([policies{j}, '_x']){n}(:, N*(h-1) + k) = zeros(J, 1);
                    continue; 
                end
                p.([policies{j}, '_x']){n}(:, N*(h-1) + k) = ...
                    (sim_xp{h}.(policies{j}){N-(k-n)}(:) - sim0.(policies{j}){N-(k-n)}(:)) / hx(h*N);
            end
        end
    end
end


%% JACOBIANS: DISTRIBUTION
ssAT = (speye(J) + dt * ss.AT);

[p.g_z, p.g_x] = deal(cell(N, 1));
for k = 1:N
    p.g_z{1}(:, k) = zeros(J, 1);
end
for k = 1:K
    p.g_x{1}(:, k) = zeros(J, 1);
end

for k = 1:N
    for n = 2:N
        p.g_z{n}(:, k) = ssAT * p.g_z{n-1}(:, k) ...
            + dt * ss.Da' * (p.s_z{n-1}(:, k) .* ss.g(:));
    end
end

for k = 1:K
    for n = 2:N
        p.g_x{n}(:, k) = ssAT * p.g_x{n-1}(:, k) ...
            + dt * ss.Da' * (p.s_x{n-1}(:, k) .* ss.g(:));
    end
end


%% JACOBIANS: MARKET CLEARING
% diff_L = X(N+1:end) - sum(sum(param.zz .* param.u1(c_dense) .* sim.g{n} .* G_dense.dx));
ss.u1z = param.zz .* param.u1(ss.c);
ss.u2z = param.zz .* param.u2(ss.c);

H.H_z = zeros(K, N);
for n = 1:N
    for k = 1:N
        H.H_z(n, k) = (ss.s(:)' * p.g_z{n}(:, k) + ss.g(:)' * p.s_z{n}(:, k)) * G_dense.dx;
        H.H_z(N+n, k) = - (ss.u1z(:)' * p.g_z{n}(:, k) + ss.g(:)' * (ss.u2z(:) .* p.c_z{n}(:, k))) * G_dense.dx;
    end
end

H.H_x = zeros(K, K);
for n = 1:N
    for k = 1:K
        H.H_x(n, k) = (ss.s(:)' * p.g_x{n}(:, k) + ss.g(:)' * p.s_x{n}(:, k)) * G_dense.dx;
        if n == k-N
            H.H_x(N+n, k) = 1 - (ss.u1z(:)' * p.g_x{n}(:, k) + ss.g(:)' * (ss.u2z(:) .* p.c_x{n}(:, k))) * G_dense.dx;
        else
            H.H_x(N+n, k) = - (ss.u1z(:)' * p.g_x{n}(:, k) + ss.g(:)' * (ss.u2z(:) .* p.c_x{n}(:, k))) * G_dense.dx;
        end
    end
end

end

