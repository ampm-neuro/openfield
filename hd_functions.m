%hd functions

%all firing rates and instantaneous HDs and pos
[all_counts_bw, all_countsz_ia, context_bw, context_ia,...
all_counts_bw_z, all_countsz_ia_z, region_context_bw, region_context_ia,...
all_pos_bw, all_hd_bw, all_pos_ia, all_visit_order_bw, all_visit_order_ia,...
all_hd_ia, all_obj_ia] = ALL_countwindow(.25);

%calculates all rate distributions and makes figures. takes about 72 hours
%to run on Dell. Auto saves variables on a loop.
%ALL_HD_relative_to_obj;

%load variables
%load('HD_rel_obj_distribs.mat')

%most important is rate_dists, a cell with all HD rate distributions
%columns for contexts Arrg1 Arrg2 Appear1 Appear2
%rows for objects 1-4 and allocentric HD

%computes correllations between rate distributions
%[withins, alls, object_corrs_alo, obj_appear] = corr_rate_dists(rd);

%plots and correllations for BW contexts
%{
[rate_dists, corrs] = ALL_HD_ctxt;
figure; hold on; 
plot(corrs(:,1), corrs(:,2), 'ko')
plot([-1 1], [-1 1], 'k--')
axis square;axis([-.5 1 -.5 1])
xlabel('Within context (r)')
ylabel('Between context (r)')
title('Similarity between HD firing rate distributions')
%}

%calculates 
[mis_all_allo, mis_all_obj, rays_all_allo, rays_all_obj, rs_all_allo, rs_all_obj] = ALL_HD_scores(rd, allo_hd);


%ALL_dist_relative_to_obj;
[all_rd_heat] = ALL_combined_distHDobj;

%calculates difference between observed and expected firing rates for each
%object given the cell's objHD tuning curve. can play a video
[expected_rates, FRs, error_FRs, ObjHD, ObjDistance, fit_objHDdistrib] = ...
    objHDtune_fit(neuron, context, rd(aidx), rs_all_obj(aidx,:), ...
    context_ia, all_counts_ia, all_pos_ia, all_hd_ia, all_obj_ia);

%plots a single frame from above video
objHDtune_fit_fig(neuron, context, frm, rd(aidx), rs_all_obj(aidx,:), context_ia, all_counts_ia, all_pos_ia, all_hd_ia, all_obj_ia);


%calculates r values corelating expected and oberseved firing rate for each
%context
rp_normmean = nan(sum(aidx),4,2); 
for a_in = 1:size(rp_normmean,1)
    for c_in = 1:4
        [~, ~, ERs_sumnorm, FRs_norm] = objHDtune_fit(a_in,...
            4+c_in, rd(aidx), rs_all_obj(aidx,:), context_ia,...
            all_counts_ia, all_pos_ia, all_hd_ia, all_obj_ia); 
        [r,p] = corr(ERs_sumnorm, FRs_norm'); rp_normmean(a_in, c_in, :) = [r,p]; 
    end
end
cell_1{1} = mean(rp_normmean(rs_all_allo(aidx,2)>0.05 & rs_all_obj(aidx,2)<=0.05, :, 1),2);
cell_1{2} = mean(rp_normmean(rs_all_allo(aidx,2)<=0.05 & rs_all_obj(aidx,2)>=0.05, :, 1),2);
cell_1{3} = mean(rp_normmean(rs_all_allo(aidx,2)>0.05 & rs_all_obj(aidx,2)>0.05, :, 1),2);


%to estimate heatmap
%go pxl by pxl, calculate tuning curve of each obj at 1:360 deg and
%average.
for neuron = [2]
    objHDtune_space(neuron, rd, rs_all_obj, rate_matrices_ia, hd_proportions_ia)
end
