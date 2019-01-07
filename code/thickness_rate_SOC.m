% len_c_ratios = linspace(0.5, 2, 7);
% len_a_ratios = linspace(1/3, 6/3, 6);
% disc_rates = logspace(-0.1, 0.5, 6);
len_c_ratios = linspace(2, 0.5, 13);
len_a_ratios = linspace(6/3, 1/3, 11);
disc_rates = logspace(-0.1, 0.5, 10) * (-1);

SOC_c_s =   zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));
SOC_a_s =   zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));
outputV_s = zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));
er_s =      zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));

for i = 1:numel(len_c_ratios)
    for j = 1:numel(len_a_ratios)
        for k = 1:numel(disc_rates)
            [ SOC_c_s(i,j,k), SOC_a_s(i,j,k), outputV_s(i,j,k), er_s(i,j,k) ] = init_fun_s(len_c_ratios(i),len_a_ratios(j), disc_rates(k));
        end
    end
end
     
%%
lc =   zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));
la =   zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));
dr =   zeros(numel(len_c_ratios), numel(len_a_ratios), numel(disc_rates));

for i = 1:numel(len_c_ratios)
    for j = 1:numel(len_a_ratios)
        for k = 1:numel(disc_rates)
            lc(i,j,k) = len_c_ratios(i);
            la(i,j,k) = len_a_ratios(j);
            dr(i,j,k) = disc_rates(k);
        end
    end
end

lc = lc(:); la = la(:); dr = dr(:); 
er_s = er_s(:); SOC_c_s = SOC_c_s(:); SOC_a_s = SOC_a_s(:);

figure;
uer = [11 31 34 41]; %uer = unique(er_s);
for i = 1:numel(uer)
    idx = er_s == uer(i);
    scatter3(lc(idx), la(idx), dr(idx), 10, 'filled');
    hold all;
end

ax = gca;
ax.XDir = 'normal';
set(gca,'zscale','log');
xlabel('relative cathode thickness');
ylabel('relative anode thickness');
zlabel('rate');
title('Exit reason');
legend('electrolyte depletion',...
    'cutoff SOC at cathode',...
    'cutover SOC at anode',...
    'cutoff voltage');

figure;
scatter3(lc, la, dr, 10, SOC_c_s, 'filled');
ax = gca;
ax.XDir = 'normal';
set(gca,'zscale','log');
xlabel('relative cathode thickness');
ylabel('relative anode thickness');
zlabel('rate');
title('Cathode SOC at the end');

figure;
scatter3(lc, la, dr, 10, SOC_a_s, 'filled');
ax = gca;
ax.XDir = 'normal';
set(gca,'zscale','log');
xlabel('relative cathode thickness');
ylabel('relative anode thickness');
zlabel('rate');
title('Anode SOC at the end');

%% Plot the results

% h = figure;
% % plot(out1.time{1},out1.Voltage{1},'--','LineWidth',2)
% A = axes;
% mesh(disc_rates, len_p_ratios, SOC_c_s);
% set(A,'XScale','log');
% ylabel('relative length of positive electrode (LTO)'); 
% xlabel('discharge rates');
% zlabel('remained SOC (graphite)');
% grid on
% box on