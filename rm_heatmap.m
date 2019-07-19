function rm_fig = rm_heatmap(rm)
%plotting funciton for distance from obj X HD rel to obj rate matrices 
% output by ALL_combined_distHDobj.m
    rm_fig = rm(:, :);
    %rm_fig = rm_fig(:,[floor(size(rm_fig,2)/2)+1:size(rm_fig,2) 1:floor(size(rm_fig,2)/2)]);
    rm_fig = smooth2a(rm_fig,1);
    
    for i = 1:size(rm_fig,1)
        rm_fig(i,:) = zscore_mtx(rm_fig(i,:)')'; 
    end
    
    rm_fig = csmooth(rm_fig);
    
    figure; imagesc(rm_fig);
    axis square; set(gca,'TickLength',[0, 0]); box off
    ylabel('Distance from object (m)')
    xlabel('Head direction (degrees from object)')
    xticks([1 size(rm_fig,2)/2 size(rm_fig,2)])
    xticklabels({'-180' 'Object' '180'});
    yticks([1 size(rm_fig,1)])
    yticklabels({'1' '0'});
    colorbar
    caxis([-1.5 1.5])
end


function mtx = csmooth(mtx)
    pages = nan(size(mtx,1), size(mtx,2), size(mtx,2));
    nanidx = isnan(mtx);
    for i = 1:size(pages,2)
        hold_ = mtx(:,[i:end 1:i-1]);
        hold_(nanidx) = nan;
        hold_ = smooth2a(hold_,3);
        pages(:,:,i) = hold_(:,[end-(i-2):end 1:end-(i-1)]);
    end
    mtx = nanmean(pages,3);
end