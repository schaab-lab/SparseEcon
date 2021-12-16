function [basis_sim, basis_agg] = estimation_model(sim, agg, type, agg_dims)

switch type   
    case 'linear'
        basis_sim = [ones(sim.N, 1), sim.X];
        basis_agg = [ones(agg.J, 1), agg.X(:, agg_dims)];

    case 'quadratic'
        basis_sim = [ones(sim.N, 1), sim.X, sim.X.^2];
        basis_agg = [ones(agg.J, 1), agg.X(:, agg_dims), agg.X(:, agg_dims).^2];
        
    case 'quadratic1'
        basis_sim = [ones(sim.N, 1), sim.X, sim.X(:, 1).^2];%, prod(sim.X, 2)];
        basis_agg = [ones(agg.J, 1), agg.X(:, agg_dims), agg.X(:, 1).^2];%, prod(agg.X, 2)];
        
    case 'quadratic2'
        basis_sim = [ones(sim.N, 1), sim.X, sim.X(:, 2).^2];%, prod(sim.X, 2)];
        basis_agg = [ones(agg.J, 1), agg.X(:, agg_dims), agg.X(:, 2).^2];%, prod(agg.X, 2)];

    case 'quadraticCov'
        basis_sim = [ones(sim.N, 1), sim.X, sim.X.^2, prod(sim.X, 2)];
        basis_agg = [ones(agg.J, 1), agg.X(:, agg_dims), agg.X(:, agg_dims).^2, prod(agg.X(:, agg_dims), 2)];        
        
end