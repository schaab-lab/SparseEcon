function hjb = HJB(V, G, param)


%% VF DERIVATIVES
num0 = 1e-8; % numerical 0 for upwind scheme

VkF = deriv_sparse(G, V, 1, 'D1F');
VkB = deriv_sparse(G, V, 1, 'D1B');

VkF = max(VkF, num0);
VkB = max(VkB, num0);

%% UPWIND
cF = param.u1inv(VkF); cF(G.grid(:, 1) == 1) = G.income(G.grid(:, 1) == 1);
cB = param.u1inv(VkB); cB(G.grid(:, 1) == 0) = G.income(G.grid(:, 1) == 0);
c0 = G.income;

sF = G.income - cF;
sB = G.income - cB;

% Upwind #1:
HF = param.u(cF) + VkF .* sF;
HB = param.u(cB) + VkB .* sB;

IF = (sF > num0) .* ((sB >= 0) + (sB <-num0).*(HF>=HB)); % not clear whether it should be (sB>=0) for BC, or (sB>-num0) for num error
IB = (sB <-num0) .* ((sF <= 0) + (sF > num0).*(HB>=HF)); % not clear whether it should be (sF<=0) for BC, or (sF<-num0) for num error
I0 = ~IF & ~IB;

% Upwind #2:
% IF = (sF > num0);        % BC takes care of this: (G.grid(:,1)<1)
% IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:,1)>0)
% I0 = ~IF & ~IB;

s = sF.*IF + sB.*IB;
c = cF.*IF + cB.*IB + c0.*I0;
u = param.u(c);

%% OUTPUT
hjb.c = c; hjb.s = s; hjb.u = u;

end

