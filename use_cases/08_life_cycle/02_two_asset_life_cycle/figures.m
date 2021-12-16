
figure; scatter3(G.a, G.t, hjb.c(:, 1));
figure; scatter3(G.a, G.t, hjb.s(:, 1));
figure; scatter3(G.a, G.t, hjb.c(:, 2));
figure; scatter3(G.a, G.t, hjb.s(:, 2));

figure; scatter3(G.k, G.t, hjb.iota(:, 1));
figure; scatter3(G.k, G.t, hjb.iota(:, 2));


figure; scatter3(G.a, G.t, g(:, 2));

figure; scatter3(G.k, G.t, g(:, 1));
figure; scatter3(G.k, G.t, g(:, 2));


figure; scatter3(G.a, G.t, V(:, 2));


figure; hold on; 
plot(G.a(G.t == param.tmax), G.V0(G.t == param.tmax, 2)); 
plot(G.a(G.t == param.tmax), 1e-8 * param.u(1e-8 + G.a(G.t == param.tmax))); hold off;


