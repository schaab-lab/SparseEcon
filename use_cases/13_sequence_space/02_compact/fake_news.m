function [H, p] = fake_news(x0, c0, z0, ss, G, G_dense, param)


%% SETUP
f = @(x, z) transition(x, z, ss, G, G_dense, param, 'fake_news');

sim0 = f(x0, z0);

hx = param.phi_jacobian * sqrt(eps) * max(abs(x0(:)), 1);
hz = param.phi_jacobian * sqrt(eps) * max(abs(z0(:)), 1);

N = param.N; K = numel(x0); J = param.discrete_types * G.J; dt = param.dt;

policies = {'s', 'c', 'm'};
aggregates = {'S', 'M'};

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


%% EQUILIBRIUM MAP
f = @(x, c, z) equilibrium_map(x, c, z, ss, aggregates, param);
obj0 = f(x0, c0, z0);

hx = param.phi_jacobian * sqrt(eps) * max(abs(x0(:)), 1);
hc = param.phi_jacobian * sqrt(eps) * max(abs(c0(:)), 1);
hz = param.phi_jacobian * sqrt(eps) * max(abs(z0(:)), 1);

% H_x:
Jx = zeros(numel(obj0), numel(x0));
parfor k = 1:numel(x0)
    x1    = x0;
    x1(k) = x1(k) + hx(k);
    obj1  = f(x1, c0, z0);

    delta  = obj1(:) - obj0(:);
    Jx(:,k) = delta / hx(k);
end

% H_c:
Jc = zeros(numel(obj0), numel(c0));
parfor k = 1:numel(c0)
    c1    = c0;
    c1(k) = c1(k) + hc(k);
    obj1  = f(x0, c1, z0);

    delta  = obj1(:) - obj0(:);
    Jc(:,k) = delta / hc(k);
end

% H_z:
Jz = zeros(numel(obj0), numel(z0));
parfor k = 1:numel(z0)
    z1    = z0;
    z1(k) = z1(k) + hz(k);
    obj1  = f(x0, c0, z1);

    delta  = obj1(:) - obj0(:);
    Jz(:,k) = delta / hz(k);
end

% Construct C_x and C_z:
C_z = zeros(numel(c0), numel(z0));
for j = 1:numel(aggregates)
    for n = 1:N
        for k = 1:numel(z0)
            % C_z(n, k) = ss.c' * p.g_z{n}(:, k) + ss.g' * p.c_z{n}(:, k);
            C_z(n+(j-1)*N, k) = (ss.([lower(aggregates{j})])(:)' * p.g_z{n}(:, k) ...
                + ss.g(:)' * p.([lower(aggregates{j}), '_z']){n}(:, k)) * G_dense.dx;
        end
    end
end

C_x = zeros(numel(c0), numel(x0));
for j = 1:numel(aggregates)
    for n = 1:N
        for k = 1:numel(x0)
            C_x(n+(j-1)*N, k) = (ss.([lower(aggregates{j})])(:)' * p.g_x{n}(:, k) ...
                + ss.g(:)' * p.([lower(aggregates{j}), '_x']){n}(:, k)) * G_dense.dx;
        end
    end
end

% Total derivative:
H.H_x = Jx + Jc * C_x;
H.H_z = Jz + Jc * C_z;


end

