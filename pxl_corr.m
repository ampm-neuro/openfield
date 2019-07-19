function [cell_pixle_corr, within_corrs] = pxl_corr(eptrials, clust, sel_contexts)
%finds correllation within and between, for example, White and Black contexts
%
%abab_order is the A B A B presentation. 
% W-B-W-B would be 1 2 3 4
% B-W-B-W would be 2 1 4 3
% A1-O1-A2-O2 would be 1 2 3 4
% B-W-O1-O2-A2-A1-B would be 6 3 5 4
%
%cell_pixle_corr is a 2 item vector. (1) is within correlation, (2) is
%between

    matrices = heatmaps_trial(eptrials, clust, sel_contexts);
    
    %poor man's reshape
    matrices_vectors = nan(numel(matrices(:,:,1)), 4);
    for i = 1:size(matrices,3)
        hold = matrices(:,:,i); hold = hold(:);
        matrices_vectors(:,i) = hold;
    end
    
    %pairwise common pixles indices
    function cp = pcpi(vct_nums, matrices_vectors)
        cp = ~isnan(matrices_vectors(:,vct_nums(1))) & ~isnan(matrices_vectors(:,vct_nums(2)));
    end
    cp_12 = pcpi([1 2], matrices_vectors);
    cp_13 = pcpi([1 3], matrices_vectors);
    cp_24 = pcpi([2 4], matrices_vectors);
    cp_34 = pcpi([3 4], matrices_vectors);
    
    %correlations

    within_corrs = [corr(matrices_vectors(cp_13,1), matrices_vectors(cp_13,3))...
        corr(matrices_vectors(cp_24,2), matrices_vectors(cp_24,4))]; %arrange, ID
    
    corr_within = mean(within_corrs);
    corr_between = mean([corr(matrices_vectors(cp_12,1), matrices_vectors(cp_12,2))...
        corr(matrices_vectors(cp_34,3), matrices_vectors(cp_34,4))]);

    %output
    cell_pixle_corr = [corr_within corr_between];

end