function hjb = HJB(V, G, param)


%% VF DERIVATIVES
num0 = 1e-10; % numerical 0 for upwind scheme

[VnF, VnB, Vnn] = deal(zeros(G.J, param.discrete_types));
for j = 1:param.discrete_types
    VnF(:,j) = deriv_sparse(G, V(:,j), 1, 'D1F', num2str(j));
    VnB(:,j) = deriv_sparse(G, V(:,j), 1, 'D1B', num2str(j));
    Vnn(:,j) = deriv_sparse(G, V(:,j), 1, 'D2',  num2str(j));
end

VnF = max(VnF, num0);
VnB = max(VnB, num0);


%% PORTFOLIO CHOICE
IP = G.n>=0; IN = G.n<0;
thetaF = - VnF ./ (G.n .* Vnn) .* (G.muR-G.r)./param.sigR2;
% constraint: theta>0
thetaF = IP.*max(thetaF,0) + IN.*min(thetaF,0);
% constraint: theta < 1 - nmin/n
thetaF = IP.*min(thetaF, 1 - param.nmin./G.n) + IN.*max(thetaF, 1 - param.nmin./G.n);    

thetaB = - VnB ./ (G.n .* Vnn) .* (G.muR-G.r)./param.sigR2;
% constraint: theta>0
thetaB = IP.*max(thetaB,0) + IN.*min(thetaB,0);
% constraint: theta < 1 - nmin/n
thetaB = IP.*min(thetaB, 1 - param.nmin./G.n) + IN.*max(thetaB, 1 - param.nmin./G.n);

% 0-case
theta0 = (thetaB+thetaF)/2; 

if any(any(isnan([thetaF, thetaB, theta0]) | isinf([thetaF, thetaB, theta0])))
    error('NaNs or Infs detected in policy function.')
end


%% CONSUMPTION & SAVINGS
cF = param.u1inv(VnF);
cB = param.u1inv(VnB);
c0 = G.income + (G.muR-G.r).*theta0.*G.n;

sF = G.income + (G.muR-G.r).*thetaF.*G.n - cF; 
sB = G.income + (G.muR-G.r).*thetaB.*G.n - cB;

% No longer restrict IF < nmax! With asset-pricing BC, household can save at nmax! No longer use Huggett BC at nmax!
IF = (sF > num0); 
IB = (sB <-num0) & ~IF;
I0 = ~IF & ~IB;


%% OUTPUT
theta = thetaF.*IF + thetaB.*IB + theta0.*I0;
s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

Vn0 = param.u1(c0);
Vn  = VnF.*IF + VnB.*IB + Vn0.*I0;

% Collect output: 
hjb.c = c; hjb.s = s; hjb.u = u; hjb.theta = theta; hjb.Vn = Vn; hjb.Vnn = Vnn;


%% CHECKSUM

% BOUNDARY CONDITION CHECK:

% Huggett BC at nmin:
assert(max(max(abs(VnB(G.n == param.nmin, :) - param.u1(G.income(G.n == param.nmin, :))))) <= num0);

% Huggett BC at nmax:
% => We NO LONGER USE the Huggett BC at nmax! We use the asset pricing BC
% instead, which *endogenously* solves for VnF!
% assert(max(max(abs(VnF(G.n==param.nmax, :) ...
%     - param.u1(G.income(G.n==param.nmax) - param.nmax * (G.muR-G.r).^2 / (param.gamma*param.sigR2)) ))) <= num0);

% Asset pricing BC at nmax:
assert(max(max(abs(Vnn(G.n == param.nmax, :) - -param.gamma / param.nmax * VnF(G.n == param.nmax, :)))) <= num0);


% POLICY FUNCTION CHECK:

% Huggett BC at nmin must still guarantee no dissaving:
assert(all(all( sB(G.n == param.nmin, :) >= -num0 )));
assert(all(all( s(G.n == param.nmin, :)  >= -num0 )));

% Budget constraint:
assert( max(max(abs( s - G.income - (G.muR-G.r).*theta.*G.n + c ))) < num0 );

% Portfolio choice: 
assert(all( (sign(theta(:,1)) == sign(G.n) | theta(:,1) == 0) ));
assert(all( (sign(theta(:,2)) == sign(G.n) | theta(:,2) == 0) ));
assert(all(all( theta(G.n == param.nmin, :) == 0 )));

% These are now implemented via asset-pricing BC at nmax:
% thetaF(G.n == param.nmax, :) = (G.muR-G.r) / (param.gamma*param.sigR2); 
% thetaB(G.n == param.nmax, :) = (G.muR-G.r) / (param.gamma*param.sigR2); 

% We no longer fix thetaB at nmax! New asset-pricing BC only applies to thetaF.
% assert(all(all( theta(G.n == param.nmax, :) == (G.muR-G.r)./(param.gamma*param.sigR2) )));

% No longer exactly 0 because BC is used, so num0 error is introduced:
assert( max(max(abs( thetaF(G.n == param.nmax, :) - (G.muR-G.r)./(param.gamma*param.sigR2) ))) < num0 );

assert(all(all( theta(G.n>0, :) <= (G.n(G.n>0) - param.nmin) ./ G.n(G.n>0) + num0 )));
assert(all(all( theta(G.n<0, :) >= (G.n(G.n<0) - param.nmin) ./ G.n(G.n<0) - num0 )));

% Leverage constraint: (this is redundant - always guaranteed given above)
assert(all(all( thetaF(G.n == param.nmin, :) == 0 )));
assert(all(all( thetaB(G.n == param.nmin, :) == 0 )));
assert(all(all( theta0(G.n == param.nmin, :) == 0 )));


%% ALTERNATIVE OUTPUT
% This produces numerically indistinguishable results!
% Vn0 = param.u1(c0);
% Vn  = VnF.*IF + VnB.*IB + Vn0.*I0;
% 
% theta2 = - Vn ./ (G.n .* Vnn) .* (G.muR-G.r)./param.sigR2;
% theta2 = IP.*max(theta2, 0) + IN.*min(theta2, 0);
% theta2 = IP.*min(theta2, 1 - param.nmin./G.n) + IN.*max(theta2, 1 - param.nmin./G.n);
% theta2(G.n == param.nmax, :) = (G.muR-G.r) / (param.gamma*param.sigR2); 
% theta2(G.n == param.nmin, :) = 0;
% 
% theta = theta2;
% c = cF.*IF + cB.*IB + c0.*I0;
% s = G.income + (G.muR-G.r).*theta.*G.n - c;
% u = param.u(c);


end


