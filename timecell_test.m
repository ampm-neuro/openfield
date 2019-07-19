%plot sorted heatmaps of timecell peaks in each context

%for each context
for ictx = 1:4
   
    %firing rates from that context
    cctx_rates = all_counts_bw(context_bw==ictx, :);
    %cctx_rates = all_countsz_ia(context_ia==ictx, :);
    %cctx_rates = cctx_rates(1:400,:);%debug
    
    cctx_rates = zscore_mtx(cctx_rates); %zscoreds
    
    %find all maximums
    all_max = nan(1,size(cctx_rates,2));
    for icell = 1:size(cctx_rates,2)
        all_max(icell) = find(cctx_rates(:,icell) == max(cctx_rates(:,icell)), 1, 'first');
    end
    
    %sorted
    [~,sort_idx] = sort(all_max);
    cctx_rates = cctx_rates(:,sort_idx);
    
    figure
    imagesc(cctx_rates')
    caxis([-3 3])
    colorbar
end