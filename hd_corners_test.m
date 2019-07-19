%edit inside hd_cell to change boundaries that bound rat to corner or
%center of arena
[rate_distributions_CT, hd_counts_CT] = hd_cell(eptrials(ismember(eptrials(:,6),5:8),:), clusters(:,1), 11, 7, 1);
[rate_distributions_SW, hd_counts_SW] = hd_cell(eptrials(ismember(eptrials(:,6),5:8),:), clusters(:,1), 11, 7, 1);
[rate_distributions_NE, hd_counts_NE] = hd_cell(eptrials(ismember(eptrials(:,6),5:8),:), clusters(:,1), 11, 7, 1);

figure; hold on
plot(zscore_mtx(rate_distributions_CT(2,:)), '-', 'Color', [0    0.4470    0.7410])
for i = 0:90:270; plot([45+i 45+i], ylim, '-', 'Color', [0    0.4470    0.7410]); end %box corners

plot(zscore_mtx(rate_distributions_SW(2,:)), '-', 'Color', [0.8500    0.3250    0.0980])
plot([45+120 45+120], ylim, '-', 'Color', [0.8500    0.3250    0.0980])%a box corner from this vantage
plot([45+240 45+240], ylim, '-', 'Color', [0.8500    0.3250    0.0980])%a box corner from this vantage

plot(zscore_mtx(rate_distributions_NE(2,:)), '-', 'Color', [0.9290    0.6940    0.1250])
plot([45+60 45+60], ylim, '-', 'Color', [0.9290    0.6940    0.1250])%a box corner from this vantage
plot([45+60+240 45+60+240], ylim, '-', 'Color', [0.9290    0.6940    0.1250])%a box corner from this vantage