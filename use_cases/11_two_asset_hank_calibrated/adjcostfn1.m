function f = adjcostfn1(iota, k, param)

f = (  param.psi0 + param.psi1 * param.psi2 * (abs(iota) ./ max(k, param.psi3)) .^ (param.psi2 - 1)) .* (iota > 0) + ...
    ( -param.psi0 - param.psi1 * param.psi2 * (abs(iota) ./ max(k, param.psi3)) .^ (param.psi2 - 1)) .* (iota < 0);

end