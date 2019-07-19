idx = [[rs_all_obj(:,2)>0.05 & rs_all_allo(:,2)>0.05]';
    [rs_all_obj(:,2)<0.05]';
    [rs_all_obj(:,2)>0.05 & rs_all_allo(:,2)<0.05]';
    [rs_all_obj(:,2)>0.05]'];


for i1 = 1:4

    pna = peaks_new_allobj(idx(i1,:), :);
    mean_heading = nan(size(pna,1),1);
    for i= 1:size(pna,1)
        %mean_heading(i) = rad2deg(circ_mean(deg2rad(pna(i,:)')));
        mean_heading(i) = pna(i,1);
    end
    mean_heading(mean_heading<0) = mean_heading(mean_heading<0) + 360;
    [~,srt_idx] = sort(mean_heading);
    pna = pna(srt_idx,:);
        
    
    figure
    for i2 = 1:4
        subplot(1,4, i2)
        imagesc(pna(:, [1 2 3 4] + 4*(i2-1))); colormap hsv
        set(gca,'TickLength',[0, 0]); box off
        axis off
    end
end