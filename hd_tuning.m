function [peaks, peaks_idx] = hd_tuning(rd)
%plots hd tuning of each cell

colors =       [0    0.4470    0.7410; %blue arrg 1
    0.4660    0.6740    0.1880; %grn arrg 2
    0.6350    0.0780    0.1840; %red obj 1
    0.9290    0.6940    0.1250]; %yellow obj 2

peaks = nan(length(rd),1); 
peaks_idx = nan(size(peaks));
peaks_local = nan(16,1); 
peaks_local_idx = nan(size(peaks_local));

for icell = 1:161 
    count = 0; 
    figure; 
    for icol = 1:4 
        for irow = 1:4 
            count = count+1; 
            
            %rotate and zscore FR distribution
            rot_idx = [size(rd{icell}{irow,icol},2)/2+1:size(rd{icell}{irow,icol},2)...
                1:size(rd{icell}{irow,icol},2)/2];
            hd_dist_hold(count,:) = rd{icell}{irow,icol}(rot_idx); 
            hd_dist_hold(count,:) = zscore_mtx(hd_dist_hold(count,:)'); 
            
            %plot relative to obj
            subplot(1, 2, 2)
            
            c_rd = hd_dist_hold(count,:);  
            if ismember(icol,2:4) || ismember(irow,2:4) 
                hold on
            end
            
            polarplot(deg2rad(1:3:360), c_rd-repmat(min(c_rd), size(c_rd)), 'color', colors(icol,:))
            ax=gca;
            ax.ThetaZeroLocation = 'top';
            ax.ThetaDir = 'clockwise';
            ax.ThetaTick = [0 90 180 270];
            ax.RTick = [];
            title('Relative to Object')
            
            %ylim
            if irow==1 && icol==1
                ax.RLim = [0 max(c_rd)-min(c_rd)];
            elseif max(ax.RLim) < max(c_rd)-min(c_rd)
                ax.RLim = [0 max(c_rd)-min(c_rd)];
            end

            %{
            %line plot
            hold on
            plot(hd_dist_hold(count,:), 'color', colors(icol,:)); 
            xlim([1 120]); 
            xticks([1 30:30:120]); 
            xticklabels({'-179' '-90' 'Obj' '90' '180'}); 
            set(gca,'TickLength',[0, 0]); box off
            %}
            
            %load peak local
            max_idx = find(hd_dist_hold(count,:)==max(hd_dist_hold(count,:)));
            peaks_local(count) = rot_idx(max_idx);
            peaks_local_idx(count) = max_idx;

        end

        rot_idx = [size(rd{icell}{5,icol},2)/2+1:size(rd{icell}{5,icol},2)...
        1:size(rd{icell}{5,icol},2)/2];
        hd_dist_hold_alo(count,:) = rd{icell}{5,icol}(rot_idx); 
        hd_dist_hold_alo(count,:) = zscore_mtx(hd_dist_hold_alo(count,:)');
        c_rd_allo = hd_dist_hold_alo(count,:); 
        
        %plot allocentric hd
        subplot(1, 2, 1)
        if ismember(icol,2:4)
            hold on
        end
        polarplot(deg2rad(1:360), c_rd_allo-repmat(min(c_rd_allo), size(c_rd_allo)), 'color', colors(icol,:))
        ax=gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        ax.ThetaTick = [0 90 180 270];
        ax.RTick = [];
        title('Allocentric')
        


        %ylim
        if icol==1
            ax.RLim = [0 max(c_rd_allo)-min(c_rd_allo)];
        elseif max(ax.RLim) < max(c_rd_allo)-min(c_rd_allo)
            ax.RLim = [0 max(c_rd_allo)-min(c_rd_allo)];
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
    peaks(icell) = circ_mean(deg2rad(peaks_local));
    peaks_idx(icell) = circ_mean(deg2rad(peaks_local_idx))
    
    %save
    %{
    var_name = ['HD_rel_obj_' num2str(icell)];
    print(['C:\Users\ampm1\Desktop\context paper\HD_rel_obj\overlays_polar\' var_name],...
        '-dpdf', '-painters', '-bestfit')
    %}
    close
    
    
end
peaks = rad2deg(peaks);
peaks_idx = rad2deg(peaks_idx);

%
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
