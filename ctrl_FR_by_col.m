function epout = ctrl_FR_by_col(eptrials, clusters, col, bins, tw_size, slide_size)
%output FRs that have been processed to remove effect of the variable in
%col

%   eptrials(:,7) contains the rat's velocity.
%   eptrials(:,8) contains the rat's acceleration
%   eptrials(:,9) contains the rat's head direction.

    %sliding time window
    %tw_size = 500; %ms
    %slide_size = 500; %ms

%for speed
lclusters = length(clusters);

%MODEL EFFECT OF COL
%
    %calculate binned rates
    rates_all = [];
    col_all = [];
    means_all = nan(4, size(clusters,1));
    stds_all = nan(4, size(clusters,1));
    for ctxt = 1:4

        %bin FRs by col
        cur_eptrials = eptrials(eptrials(:,5)==ctxt, :);
        [rates_out, col_out] = FR_by_col(cur_eptrials, clusters, col, bins, tw_size, slide_size);

        %save means and stds for later
        means_all(ctxt, :) = nanmean(rates_out, 1);
        stds_all(ctxt, :) = nanstd(rates_out, 1);

        %zscore each context independently
        rates_out = zscore_mtx(rates_out);

        %combine normed contexts
        rates_all = [rates_all; rates_out];
        col_all = [col_all; col_out];
    end
    col_all = col_all;


    %model each cell
    if col == 7
        degree = 2;
    elseif col == 8
        degree = 4;
    end
    poly = nan(lclusters, 1+degree);
    FR_target = nan(lclusters,1);
    for ic = 1:lclusters
    
        %model relationship between zFR and col
        nnan_idx = ~isnan(rates_all(:,1));
        poly(ic,:) = polyfit(col_all(nnan_idx), rates_all(nnan_idx,ic), degree);

        %target FR for col variable
        if col == 7
            %if velocity, recalibrate to v = 0.35 m/s
            FR_target(ic) = polyval(poly(ic,:), 0.35);
        elseif col == 8
            %if velocity, recalibrate to a = 0.4 m/s^2
            FR_target(ic) = polyval(poly(ic,:), 0.4);
        end
    end

%RECALCULATE FRs at every TW
%
epout = eptrials;
epout_cols = size(epout, 2); %for speed
for ctxt = 1:4

    ctxt
    
    %constrain eptrials
    cur_epout = epout(epout(:,5)==ctxt, :);
    cur_ctxt = mode(cur_epout(:, 6));

    %all times
    all_times = cur_epout(cur_epout(:,4)==1,1);

    %preindex
    vid_idx = cur_epout(:,4)==1;
                
    for tw = 1:(slide_size/10):length(cur_epout(cur_epout(:,4)==1,1))-(tw_size/10)

        %tw index
        tw_idx = cur_epout(:,1)>all_times(tw) & cur_epout(:,1)<=all_times(tw+(tw_size/10));

        %current col_val
        col_val = mean(cur_epout(tw_idx & vid_idx, col));
        
        %col correction (how FR needs to be adjusted for this col val)
        current_polyval = nan(lclusters,1);
        for ic = 1:lclusters
            current_polyval(ic) = polyval(poly(ic,:), col_val);
        end
        col_correction = FR_target - current_polyval;

        %count spikes in window (sets upper bin past all clusters)
        dummy_clusters = [clusters; clusters(end)+.0001];
        CSC = histcounts(cur_epout(tw_idx & ~vid_idx,4), dummy_clusters);
        
        %convert to FR
        CFR = CSC./(tw_size/1000);
        
        %zscore current FRs
        zCFR = (CFR - means_all(ctxt, :))./stds_all(ctxt, :);
       
        %add col correction
        zCFR = zCFR + col_correction';
        
        %new firing rates (un zscore)
        NFR = (zCFR.*stds_all(ctxt, :)) + means_all(ctxt, :);
        
        %convert back to spike counts
        NSC = round(NFR.*(tw_size/1000));
        
        %fix and prepare indices for NSC
        NSC(NSC<0) = 0; %if negative corrected FR
        cumsumNSC = cumsum(NSC);
        startNSC = cumsumNSC - NSC + ones(size(NSC));
        
        %erase old spikes
        eptwidx = epout(:,1)>all_times(tw) & epout(:,1)<=all_times(tw+(tw_size/10));
        epout(eptwidx & ismember(epout(:,4), clusters), :) = [];
        
        %prep
        if sum(NSC) < 1
            continue
        elseif isnan(sum(NSC))
            continue
        elseif isinf(sum(NSC))
            continue
        end
        NSR = nan(sum(NSC), epout_cols);
        
        %new spike times (linspaced)
        for ic = 1:lclusters
            
            %linspace spikes throughout time window
            newspike_times = linspace(all_times(tw)+0.000001, all_times(tw+(tw_size/10)), NSC(ic))';
            
            %add times and cluster ids to preallocated rows
            NSR(startNSC(ic):cumsumNSC(ic), 1) = newspike_times;
            NSR(startNSC(ic):cumsumNSC(ic), 4) = clusters(ic);
        end
        
        %add context data to NSRs
        NSR(:,5:6) = repmat([ctxt cur_ctxt], size(NSR,1), 1);
        
        %slip new rows into epout by time
        epout = sortrows([epout; NSR], 1);
        
    end
                
end 

