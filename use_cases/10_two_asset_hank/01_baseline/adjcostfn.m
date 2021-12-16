function f = adjcostfn(iota, k, param)

k = max(k,param.psi3);
x = abs(iota)./k;
f = ( param.psi0 * x + param.psi1 * x .^ param.psi2 ) .* (iota>0) ...
   +( param.psi0 * x + param.psi1 * x .^ param.psi2 ) .* (iota<0);

f = f.*k;

end