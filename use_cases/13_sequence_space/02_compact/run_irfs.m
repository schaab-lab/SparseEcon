function run_irfs(sim, ss, param)

ll = [];

figure;

% Output:
subplot(2, 3, 1); hold on;
for j = 1:numel(sim)
    ll(j) = plot(param.t, 100*(sim{j}.Y - ss.Y)/ss.Y);
end
hold off; title('Output ($\%$dev)', 'Interpreter', 'Latex', 'FontSize', 15);
xlim([0, min(15, param.T)]);

% Labor:
subplot(2, 3, 2); hold on;
for j = 1:numel(sim)
    plot(param.t, 100*(sim{j}.N - ss.N)/ss.N);
end
hold off; title('Labor ($\%$dev)', 'Interpreter', 'Latex', 'FontSize', 15);
xlim([0, min(15, param.T)]);

% Wages:
subplot(2, 3, 3); hold on;
for j = 1:numel(sim)
    plot(param.t, 100 * (sim{j}.w - ss.w) / ss.w);
end
hold off; title('Real wage ($\%$dev)', 'Interpreter', 'Latex', 'FontSize', 15);
xlim([0, min(15, param.T)]);

% Inflation:
subplot(2, 3, 4); hold on;
for j = 1:numel(sim)
    plot(param.t, 100 * sim{j}.pi);
end
hold off; title('Inflation ($\%$)', 'Interpreter', 'Latex', 'FontSize', 15);
xlabel('Time (quarters)', 'FontSize', 13, 'Interpreter', 'Latex'); xlim([0, min(15, param.T)]);

% Nominal interest rate:
subplot(2, 3, 5); hold on;
for j = 1:numel(sim)
    plot(param.t, 100 * sim{j}.i);
end
hold off; title('Interest rate $i$ ($\%$)', 'Interpreter', 'Latex', 'FontSize', 15);
xlabel('Time (quarters)', 'FontSize', 13, 'Interpreter', 'Latex'); xlim([0, min(15, param.T)]);

% Shock:
subplot(2, 3, 6); hold on;
switch param.shock_type
    case 'TFP'
        plot(param.t, 100*(sim{1}.Z - ss.Z)/ss.Z);
        title_string = 'TFP Shock';
    case 'demand'
        plot(param.t, 100*(sim{1}.rho - param.rho)/param.rho);
        title_string = 'Demand Shock';
    case 'cost-push'
        plot(param.t, 100*(sim{1}.epsilon - param.epsilon)/param.epsilon);
        title_string = 'Cost-Push Shock';
    case 'monetary'
        plot(param.t, 100*sim{1}.theta);
        title_string = 'Monetary Shock';
end
hold off;
title(title_string, 'FontSize', 15, 'Interpreter', 'Latex');
xlabel('Time (quarters)', 'FontSize', 13, 'Interpreter', 'Latex'); xlim([0, min(15, param.T)]);

set(gcf,'Position',[500,500,800,400]);
set(gcf,'renderer','Painters');

end