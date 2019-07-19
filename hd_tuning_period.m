function [peak_idx] = hd_tuning_period(rd)
%plots hd tuning of each cell

peak_idx = nan(length(rd),1);
figure; hold on

for icell = 130:161%111:161 
    count = 0; 
    for icol = 1:4 
        count = count+1;
        
        %plot allocentric hd
        rot_idx = [size(rd{icell}{5,icol},2)/2+1:size(rd{icell}{5,icol},2)...
            1:size(rd{icell}{5,icol},2)/2];
        hd_dist_hold_alo(count,:) = rd{icell}{5,icol}(rot_idx); 
        hd_dist_hold_alo(count,:) = zscore_mtx(hd_dist_hold_alo(count,:)');

    end
    
    hdh = zscore_mtx(nanmean(hd_dist_hold_alo)');
    peak_idx(icell) = find(hdh==max(hdh), 1, 'first');
    rotation_fix = 45-peak_idx(icell);
    if rotation_fix < 0 
        rot_fix = [(1-rotation_fix):length(hdh) 1:(1-rotation_fix-1)];
    else
        rot_fix = [(1-rotation_fix):length(hdh)-rotation_fix 1:(1-rotation_fix-1)];
        rot_fix(rot_fix<1) = rot_fix(rot_fix<1) + 360;
    end
    plot(hdh(rot_fix)); 

end

xlim([1 360]); 
xticks([1 90:90:360]); 
set(gca,'TickLength',[0, 0]); box off
hold on; plot([45 45], ylim, 'k-')
end
