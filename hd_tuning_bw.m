function [peaks, peaks_idx] = hd_tuning_bw(rd)
%plots hd tuning of each cell

colors = ...
    [0    0.4470    0.7410; %blue arrg 1
    0.6350    0.0780    0.1840]; %red obj 1

peaks = nan(length(rd),1); 
peaks_idx = nan(size(peaks));
peaks_local = nan(16,1); 
peaks_local_idx = nan(size(peaks_local));

for icell = 1:length(rd) 
    count = 0; 
    figure; 
    for icol = 1:4 
        count = count+1;
        rot_idx = [size(rd{icell}{1,icol},2)/2+1:size(rd{icell}{1,icol},2)...
        1:size(rd{icell}{1,icol},2)/2];
        hd_dist_hold_alo(count,:) = rd{icell}{1,icol}(rot_idx); 
        hd_dist_hold_alo(count,:) = zscore_mtx(hd_dist_hold_alo(count,:)');
        c_rd_allo = smooth(hd_dist_hold_alo(count,:),3); 
        
        if ismember(icol,1:2)
            color_idx = 1;
        else 
            color_idx = 2;
        end
        
        %plot allocentric hd
        if ismember(icol,2:4)
            hold on
        end
        
        polarplot(deg2rad(1:360), c_rd_allo-repmat(min(c_rd_allo), size(c_rd_allo)), 'color', colors(color_idx,:))
        ax=gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        ax.ThetaTick = [0 90 180 270];
        ax.RTick = [];
        title('Allocentric')

        %ylim
        if icol==1
            nanmax(c_rd_allo)
            nanmin(c_rd_allo)
            ax.RLim = [0 nanmax(c_rd_allo)-nanmin(c_rd_allo)]
        elseif max(ax.RLim) < nanmax(c_rd_allo)-nanmin(c_rd_allo)
            ax.RLim = [0 nanmax(c_rd_allo)-nanmin(c_rd_allo)]
        end
        
                
        %{
        rot_idx = [size(rd{icell}{5,icol},2)/2+1:size(rd{icell}{5,icol},2)...
            1:size(rd{icell}{5,icol},2)/2];
        hd_dist_hold_alo(count,:) = rd{icell}{5,icol}(rot_idx); 
        hd_dist_hold_alo(count,:) = zscore_mtx(hd_dist_hold_alo(count,:)');
        c_rd = hd_dist_hold(count,:); 
        plot(hd_dist_hold_alo(count,:), 'color', colors(icol,:)); 
        xlim([1 360]); 
        xticks([1 90:90:360]); 
        xticklabels([181 270 0 90 180]); 
        set(gca,'TickLength',[0, 0]); box off   
        %}
        
    end
    %peaks(icell) = circ_mean(peaks_local, [1 360]);
    %peaks_idx(icell) = circ_mean(peaks_local_idx, [1 360]);
    
    %save
    var_name = ['HD_bw_' num2str(icell)];
    print(['C:\Users\ampm1\Desktop\context paper\HD_bw\overlays_polar\' var_name],...
        '-dpdf', '-painters', '-bestfit')
    
    
end
%{
    figure; hold on
    for ip = 1:length(peaks_idx)
        plot([peaks_idx(ip) peaks_idx(ip)], ylim, 'k-')
    end
    xlim([1 120]); 
    xticks([1 30:30:120]); 
    xticklabels({'-179' '-90' 'Obj' '90' '180'}); 
    set(gca,'TickLength',[0, 0]); box off
%}
end
