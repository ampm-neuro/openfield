%time window length
time_window = 0.25;%s

%spatial bin count (1 side)
bins = 6; min_visits = 20;

%all_counts vars from 
%{
[all_counts_bw, all_countsz_ia, context_bw, context_ia,...
all_counts_bw_z, all_countsz_ia_z, region_idx_bw, region_idx_ia, ...
all_pos_bw, all_pos_ia] = ALL_countwindow(time_window);
%}

%
idx_bw = context_bw;
mtx = all_counts_bw(~ismember(idx_bw, 9:11), :);
mtx = zscore_mtx(mtx);
ctx_idx_bw = context_bw(~ismember(idx_bw, 9:11));




%classify
[p_correct, assignments, class, spatial_error] = timespace_classify(mtx, ctx_idx_bw, all_pos_bw(~ismember(idx_bw, 9:11),:,1), all_pos_bw(~ismember(idx_bw, 9:11),:,2), bins, min_visits);

%[p_correct, assignments, class, spatial_error] = timespace_classify(all_counts_bw_z, context_bw, all_pos_bw(:,:,1), all_pos_bw(:,:,2), bins, min_visits);
%[p_correct, assignments, class, spatial_error] = timespace_classify(all_countsz_ia_z, context_ia, all_pos_ia(:,:,1), all_pos_ia(:,:,2), bins, min_visits);
p_correct


%scatterplot figure
%
figure; hold on;
quarter_tw = length(assignments)/4;
quarter_bin = bins^2;
for iq = 1:3
    plot([quarter_tw*iq quarter_tw*iq], [0 4*bins^2], '-', 'color', [.85 .85 .85])
    plot([0 length(assignments)], [quarter_bin*iq quarter_bin*iq], '-', 'color', [.85 .85 .85])
end
plot(assignments, '.', 'markersize', 12, 'color', [.35 .35 .35])
plot(class, 'r-')
set(gca,'TickLength',[0, 0]); box off
xlim([0 length(assignments)+1])
axis square

%}


figure; hold on;

quarter_bin = bins^2;
for iq = 1:3
    plot([quarter_bin*iq quarter_bin*iq], [0 4*quarter_bin], '-', 'color', [.85 .85 .85])
    plot([0 4*quarter_bin], [quarter_bin*iq quarter_bin*iq], '-', 'color', [.85 .85 .85])
end

histc_plot = nan(4*bins^2, 4*bins^2);
tw_rngs = reshape(1:length(assignments), min_visits, 4*bins^2);
between_errors = 0;
within_errors = 0;
for pxl = 1:4*bins^2
   tw_rng = tw_rngs(:,pxl);
   
   plot(repmat(pxl, min_visits,1), assignments(tw_rng), 'o', 'markersize', 5, 'color', [.35 .35 .35])
   
   %within and between errors
   if ismember(pxl, 1:bins^2)
       between_errors = between_errors + sum(ismember(assignments(tw_rng), 2*bins^2+1:4*bins^2));
       within_errors = within_errors + sum(ismember(assignments(tw_rng), bins^2+1:2*bins^2))/2;
       
   elseif ismember(pxl, bins^2+1:2*bins^2)
       between_errors = between_errors + sum(ismember(assignments(tw_rng), 2*bins^2+1:4*bins^2));
       within_errors = within_errors + sum(ismember(assignments(tw_rng), 1:bins^2))/2;
       
   elseif ismember(pxl, 2*bins^2+1:3*bins^2)
       between_errors = between_errors + sum(ismember(assignments(tw_rng), 1:2*bins^2));
       within_errors = within_errors + sum(ismember(assignments(tw_rng), 3*bins^2+1:4*bins^2))/2;
       
   elseif ismember(pxl, 3*bins^2+1:4*bins^2)
       between_errors = between_errors + sum(ismember(assignments(tw_rng), 1:2*bins^2));
       within_errors = within_errors + sum(ismember(assignments(tw_rng), 2*bins^2+1:3*bins^2))/2;
       
   end

   histc_plot_hold = histcounts(assignments(tw_rng), 1:4*bins^2+1)';
   histc_plot_hold = histc_plot_hold./sum(histc_plot_hold(:));
   
   histc_plot(:,pxl) = histc_plot_hold;
    
end

axis square
hold on; plot([0 bins^2*4], [0 bins^2*4], 'r-')
set(gca,'TickLength',[0, 0]); box off
xticks([1 36 36*2 36*3 36*4])
yticks([1 36 36*2 36*3 36*4])
axis([-1 bins^2*4+2 -1 bins^2*4+2])

figure
imagesc(histc_plot)

%spatial error figure
%
figure; histogram(spatial_error, 0:1:ceil(dist([1 1], [bins; bins])), 'normalization', 'probability')


%classification error figure
%
figure
bar([within_errors between_errors])


