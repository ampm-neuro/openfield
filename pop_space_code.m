%time window length
time_window = 0.25;%s
%spatial bin count (1 side)
bins = 6; min_vis = 20;

%[out_bw, out_ia, cell_file_ids_bw, cell_file_ids_ia] = ALL_pixlecorr;
%figure; plot(out_bw(:,2), out_bw(:,1), 'k.'); hold on; plot([-1 2], [-1 2], 'k--'); hold off
%[out_bw_first, out_bw_second, out_bw_all] = ALL_pop_pixlecorr(bins);
%[corms, same_diag_corr, withins_diag_corr, betweens_diag_corr] = all_pop_pixlecorr_plots(out_bw_first, out_bw_second, out_bw_all);


%calculate time window rates
%{
[all_counts_bw, all_countsz_ia, context_bw, context_ia,...
all_counts_bw_z, all_countsz_ia_z, region_idx_bw, region_idx_ia, ...
all_pos_bw, all_hd_bw, all_pos_ia] = ALL_countwindow(time_window);
%}

%{
%calculate firing rate maps for every cell, vert vectorize, and concanonate
%[out_mtx] = ALL_heatmaps_ForClassify(bins);
%out_mtx = out_mtx.*time_window; %correct for timewindow rates
%out_mtx_z = zscore_mtx(out_mtx);
%black_tms = mean(cat(3, out_mtx(1:bins^2,:), out_mtx(1+bins^2:2*bins^2,:)),3);
%white_tms = mean(cat(3, out_mtx(1+2*bins^2:3*bins^2,:), out_mtx(1+3*bins^2:4*bins^2,:)),3);
%}


%classify tws into spatial tms bins 
%
%[p_correct_withinbetween_1, assignments_withinbetween_1, spaceerror_withinbetween_1, chance_spaceerror_1] = ...
%    ctxtspace_classify_crossbw(all_counts_bw_z, context_bw, all_pos_bw(:,:,1), all_pos_bw(:,:,2),...
%    bins, min_vis);
[all_within, all_between, all_shuffle] = ...
    ctxtspace_classify_crossbw(all_counts_bw_z, context_bw, all_pos_bw(:,:,1), all_pos_bw(:,:,2),...
    bins, min_vis);

%reverse context labels for first and second visits to each context
%{
context_bw_2 = context_bw;
context_bw_2(context_bw==1) = 2;
context_bw_2(context_bw==2) = 1;
context_bw_2(context_bw==3) = 4;
context_bw_2(context_bw==4) = 3;

[p_correct_withinbetween_2, assignments_withinbetween_2, spaceerror_withinbetween_2, chance_spaceerror_2] = ...
    ctxtspace_classify_crossbw(all_counts_bw_z, context_bw_2, all_pos_bw(:,:,1), all_pos_bw(:,:,2),...
    bins, min_vis);
%}


%make heatmap of guess in a middle bin
%classcount(assignments_withinbetween, bins, bins^2/2)

%{
[p_correct, assignments, class, correct_spaceerror, incorrect_spaceerror, chance_spaceerror] = ...
    ctxtspace_classify(all_counts_bw_z, context_bw, all_pos_bw(:,:,1), all_pos_bw(:,:,2), bins);
p_correct
    %}
