function f = adjcostfn1inv(chi, k, param)

kappa = (param.psi1 * param.psi2) ^ (1/(1-param.psi2));
f = ( kappa * ( chi-param.psi0).^(1/(param.psi2-1)) ) .* (chi > param.psi0) - ...
    ( kappa * (-chi-param.psi0).^(1/(param.psi2-1)) ) .* (chi <-param.psi0);

f = f .* max(k, param.psi3);

Imax = f > param.dmax;
f = f.*(1-Imax) + param.dmax.*(Imax);

end
