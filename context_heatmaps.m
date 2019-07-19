function context_heatmaps(eptrials, context_ids, clusts)
%makes subplot heatmaps of context visits in sorted order

for cf = clusts
                
    %{
    count = count +1
    if count < 10
        continue
    end
    %}

    bins = 10;

    figure; axis off; hold on
    %ha = tight_subplot(2,2,[.05 -.3],[.005 .005],[.005 .005]);


    vis_count = 0;
    for visit = sort(context_ids)'
        vis_count = vis_count+1;
        %current visit
        visit_idx = eptrials(:,6)==visit;

        %plot heatmap
        subplot(2,2, vis_count);
        %axes(ha(vis_count))
        [rm,~,~,~] = rate_mtx(eptrials(visit_idx,:), cf, bins); %skaggsmoothed rm
        rm = inpaint_nans(smooth2a(rm,1));
        imagesc(rm); 

        if ismember(visit, 5:8)
            hold on
            object_grid(visit, bins)
        end

        if vis_count == 1
            rm_min = min(rm(:));
            rm_max = max(rm(:));
        else     
            rm_min = min([rm_min min(rm(:))]);
            rm_max = max([rm_max max(rm(:))]);
        end

        axis square; axis off
        title(num2str(mode(eptrials(visit_idx,6))))


        if vis_count == 4
            for viz = 1:4
                %axes(ha(viz))
                subplot(2,2, viz)
                colorbar
                caxis([rm_min*1.1      rm_max*.9])
            end
        end
        clearvars cc_hi_sm cc_lo_sm

    end

    %set(gcf, 'Position', floor(get(0,'Screensize')./2.9)); % Maximize figure
end