function sim = macro_block_post(x, c, z, sim, ss, param, aggregates)


%% INITIALIZE AGGREGATES
for j = 1:numel(aggregates)
    sim.([aggregates{j}]) = reshape(c(1+(j-1)*param.N : j*param.N), [param.N, 1]);
end


%% MACRO BLOCK: POST


end

