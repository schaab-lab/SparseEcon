function f = adjcostfn1(iota, k, param)

f = (  param.psi0 + param.psi1 * param.psi2 * (abs(iota) ./ max(k,param.psi3)) .^ (param.psi2 - 1)) .* (iota>0) + ...
    ( -param.psi0 - param.psi1 * param.psi2 * (abs(iota) ./ max(k,param.psi3)) .^ (param.psi2 - 1)) .* (iota<0);
%{
f = (  param.kappa0_d + (abs(iota) ./ max(k,param.kappa3) /param.kappa1_d).^param.kappa2_d) .* (iota>0) + ...
    ( -param.kappa0_w - (abs(iota) ./ max(k,param.kappa3) /param.kappa1_w).^param.kappa2_w) .* (iota<0);
%}
end