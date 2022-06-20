function mm = moments(ss, G, G_dense, param)

%% AGGREGATES
mm.KY = ss.Q * ss.K / sum(ss.Y); % @Stacy: let's leave Q as real cap price here
% mm.AY = param.Ag_share;
% mm.borrower = sum(sum((G_dense.a<0) .* ss.g * G_dense.dx));
% mm.poorhtm = sum(sum((abs(G_dense.a)<=1e-9 & abs(G_dense.k)<=1e-9) .* ss.g * G_dense.dx));
% mm.richhtm = sum(sum((abs(G_dense.a)<=1e-9 & G_dense.k>1e-9) .* ss.g * G_dense.dx));
% mm.borrow = ss.rborrow;
% mm.GY = param.gov_share;


%% MPC & MPX
% mpc = compute_mpc(1, 1, ss.c, ss.A, G, param);
% mm.average_mpc = sum(sum( (G.BH_dense * mpc .* ss.g .* G_dense.dx)));
% mm.mpx = 3*mm.average_mpc;


% %% LIQUID
% ss.gc = (ss.g(:,1) + ss.g(:,2)) * G_dense.dx; % J*1 distribution combined 
% % PDF & CDF:
% num_a = 2^param.l_dense(1)+1;
% liquid = G_dense.a(1:num_a);
% pdf_liquid = zeros(num_a,1);
% for i=1:num_a
%     pdf_liquid(i) = sum(ss.gc .* (G_dense.a==liquid(i)));
% end
% cdf_liquid = zeros(num_a,1);
% cdf_liquid(1) = pdf_liquid(1);
% for i=1:num_a-1
%     cdf_liquid(i+1) = cdf_liquid(i) + pdf_liquid(i+1);
% end
% share_liquid = zeros(num_a,1);
% share_liquid(1) = pdf_liquid(1)*liquid(1)/ss.B*100;
% for i=1:num_a-1
% share_liquid(i+1) = pdf_liquid(i+1)*liquid(i+1)/ss.B*100 + share_liquid(i);
% end
% 
% % Percentile:
% mm.liquid_t01p = 100 - interp1(cdf_liquid, share_liquid, (1-0.001), 'spline');
% mm.liquid_t1p = 100 - interp1(cdf_liquid, share_liquid, (1-0.01), 'spline');
% mm.liquid_t10p = 100 - interp1(cdf_liquid, share_liquid, (1-0.1), 'spline');
% mm.liquid_b50p = interp1(cdf_liquid, share_liquid, 0.5, 'spline');
% mm.liquid_b25p = interp1(cdf_liquid, share_liquid, 0.25, 'spline');

% % Gini:
% dis_liquid = [share_liquid,cdf_liquid] .* (share_liquid>=0);
% share_liquid_clean = nonzeros(dis_liquid(:,1))/100;
% cdf_liquid_clean = nonzeros(dis_liquid(:,2));
% area_liquid = zeros(numel(share_liquid_clean),1);
% area_liquid(1) = 0.5*cdf_liquid_clean(1)*share_liquid_clean(1);
% for i=2:numel(share_liquid_clean)
%     area_liquid(i) = area_liquid(i-1) + (share_liquid_clean(i-1)+share_liquid_clean(i))*(cdf_liquid_clean(i)-cdf_liquid_clean(i-1))/2;
% end
% areaB_liquid = area_liquid(end);
% mm.liquid_gini = (0.5-areaB_liquid)/0.5;


% %% IllIQUID
% % PDF & CDF: 
% num_k = 2^param.l_dense(2)+1;
% illiquid = unique(G_dense.k);
% pdf_illiquid = zeros(num_k,1);
% for i=1:num_k
%     pdf_illiquid(i) = sum(ss.gc .* (G_dense.k==illiquid(i)));
% end
% cdf_illiquid = zeros(num_k,1);
% cdf_illiquid(1) = pdf_illiquid(1);
% for i=1:num_k-1
%     cdf_illiquid(i+1) = cdf_illiquid(i) + pdf_illiquid(i+1);
% end
% share_illiquid = zeros(num_k,1);
% share_illiquid(1) = pdf_illiquid(1)*illiquid(1)/ss.K*100;
% for i=1:num_k-1
% share_illiquid(i+1) = pdf_illiquid(i+1)*illiquid(i+1)/ss.K*100 + share_illiquid(i);
% end
% 
% % Percentile:
% mm.illiquid_t01p = max(100 - interp1(unique(cdf_illiquid), unique(share_illiquid), (1-0.001), 'spline'),0);
% mm.illiquid_t1p = max(100 - interp1(unique(cdf_illiquid), unique(share_illiquid), (1-0.01), 'spline'),0);
% mm.illiquid_t10p = max(100 - interp1(unique(cdf_illiquid), unique(share_illiquid), (1-0.1), 'spline'),0);
% mm.illiquid_b50p = max(interp1(unique(cdf_illiquid), unique(share_illiquid), 0.5, 'spline'),0);
% mm.illiquid_b25p = max(interp1(unique(cdf_illiquid), unique(share_illiquid), 0.25, 'spline'),0);
% % 
% % % Gini:
% dis_illiquid = [share_illiquid,cdf_illiquid] .* (share_illiquid>=0);
% share_illiquid_clean = nonzeros(dis_illiquid(:,1))/100;
% cdf_illiquid_clean = nonzeros(dis_illiquid(:,2));
% area_illiquid = zeros(numel(share_illiquid_clean),1);
% area_illiquid(1) = 0.5*cdf_illiquid_clean(1)*share_illiquid_clean(1);
% for i=2:numel(share_illiquid_clean)
%     area_illiquid(i) = area_illiquid(i-1) + (share_illiquid_clean(i-1)+share_illiquid_clean(i))*(cdf_illiquid_clean(i)-cdf_illiquid_clean(i-1))/2;
% end
% areaB_illiquid = area_illiquid(end);
% mm.illiquid_gini = (0.5-areaB_illiquid)/0.5;


%% OUTPUT
% Targeted:
fprintf('Targeted moments: \n');
fprintf('KMV/Bayer illiquid assets to GDP 2.92/2.86: %.2f \n', mm.KY);
% fprintf('KMV liquid assets to GDP 0.26: %.2f \n', mm.AY);
% fprintf('ARS gov spending to output 0.16: %.2f \n', mm.GY);
% fprintf('KMV/Bayer fraction borrowers 0.15/0.16: %.2f \n', mm.borrower);
% fprintf('KMV frac with a=0 & b=0 0.1: %.2f \n', mm.poorhtm);
% fprintf('KMV frac with a=0 & b>0 0.2: %.2f \n', mm.richhtm);
% fprintf('KV quarterly MPC 0.15: %.2f \n', mm.average_mpc);
% fprintf('KMV borrowing rate 0.08: %.2f \n', mm.borrow);

% % Not targeted:
% fprintf('Untargeted moments: \n');
% fprintf('KMV capital to GDP 2.92: %.2f \n', mm.KY);
% fprintf('KMV liquid top 0.1p 17/2.3: %.2f \n', mm.liquid_t01p);
% fprintf('KMV liquid top 1p 47/18: %.2f \n', mm.liquid_t1p);
% fprintf('KMV liquid top 10p 86/75: %.2f \n', mm.liquid_t10p);
% fprintf('KMV liquid bottom 50p -4/-3: %.2f \n', mm.liquid_b50p);
% fprintf('KMV liquid bottom 25p -5/-3: %.2f \n', mm.liquid_b25p);
% fprintf('KMV illiquid top 0.1p 12/7: %.2f \n', mm.illiquid_t01p);
% fprintf('KMV illiquid top 1p 33/40: %.2f \n', mm.illiquid_t1p);
% fprintf('KMV illiquid top 10p 70/88: %.2f \n', mm.illiquid_t10p);
% fprintf('KMV illiquid bottom 50p 3/0.1: %.2f \n', mm.illiquid_b50p);
% fprintf('KMV illiquid bottom 25p 0/0: %.2f \n', mm.illiquid_b25p);
% fprintf('KMV liquid gini 0.98/0.86: %.2f \n', mm.liquid_gini);
% fprintf('KMV illiquid gini 0.81/0.82: %.2f \n', mm.illiquid_gini);

% fprintf('ARL capital output ratio 2.23\n');
% fprintf('ARL total household wealth to GDP 3.82\n');
% fprintf('ARL gov debt value to GDP 0.46 \n');
% fprintf('ARL liquidity to GDP 0.23 \n');
% fprintf('ARL average MPC 0.55: TBC \n');
% 
% fprintf('Bayer mean liquidity to mean illiquid 0.09: %.2f \n', mm.BKH);
% fprintf('Bayer second quintile liquidity to illiquid ratio 0.33 \n');
% fprintf('Bayer Gini total wealth 0.78 \n');


end