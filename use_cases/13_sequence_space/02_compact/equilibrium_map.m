function diff = equilibrium_map(x, c, z, ss, aggregates, param)

%% MACRO BLOCK: PRE
sim = macro_block_pre(x, z, ss, param);


%% MACRO BLOCK: POST
sim = macro_block_post(x, c, z, sim, ss, param, aggregates);


%% OUTPUT
% excess_bonds = sim.B;
% excess_goods = sim.Y - sim.C;
excess_union = sim.MX - sim.M;
excess_saving = sim.S;

diff = [excess_saving; excess_union];

end

