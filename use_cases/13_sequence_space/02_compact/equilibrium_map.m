function [sim, G, G_dense] = equilibrium_map(x, c, z, ss, G, G_dense, param, query)


%% MACRO BLOCK: PRE
sim = macro_block_pre(x, z, ss, param);


%% MACRO BLOCK: POST
aggregates = {'B', 'C', 'S', 'M'};
sim = macro_block_post(x, c, z, sim, ss, param, aggregates);


%% OUTPUT
sim.excess_bonds = sim.B;
sim.excess_goods = sim.Y - sim.C;
sim.excess_union = sim.MX - sim.M;
sim.excess_saving = sim.S;

sim.diff_markets = [sim.excess_goods, sim.excess_bonds, sim.excess_union, sim.excess_saving];

if nargin > 6 && ~any(ismember(query, {'all'}))
    sim = [sim.excess_saving; sim.excess_union];
end


end

