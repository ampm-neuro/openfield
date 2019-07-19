function [spike_counts, context_id, xpos, ypos, hd, obj_vars] = slide_rate_window_ia(eptrials, clusters, window_duration)
% moves through session fiding spike counts for each window of
% window_duration (no overlap) and recording the context from which it came
%
% window_duration is in seconds
%
% intended to be called iteratively, once for each unique context, although
% could be called all at once
%
% this function was written to do a context discrimination analysis
% 

%find times surrounding each bin
all_time_bins = min(eptrials(:,1)):window_duration:max(eptrials(:,1));
time_bins_low = all_time_bins(1:(end-1));
time_bins_hi = all_time_bins(2:end);

%preallocate
spike_counts = nan(length(time_bins_low), size(clusters,1));
context_id = nan(length(time_bins_low), 1);
xpos = nan(length(time_bins_low), 1);
ypos = nan(length(time_bins_low), 1);
hd = nan(length(time_bins_low), 1);
obj_vars = nan(length(time_bins_low), 8); %4xProximity, 4xObjHD

%timesave
vid_idx = eptrials(:,4)==1;

%iterate through time windows and count spikes
for wi = 1:length(time_bins_low)
    
    time_low = time_bins_low(wi);
    time_hi = time_bins_hi(wi);
    time_idx = eptrials(:,1)>=time_low & eptrials(:,1)<time_hi;
    
    %load spike_counts
    spike_counts(wi,:) = histcounts(eptrials(time_idx,4), [clusters(:,1); clusters(end,1)+.01]);
    
    %pos
    xpos(wi) = mean(eptrials(time_idx, 2));
    ypos(wi) = mean(eptrials(time_idx, 3));
    
    if isnan(sum(eptrials(time_idx & vid_idx, 9)))
        hd(wi) = nan;
    else
        hd(wi) = rad2deg(circ_mean(deg2rad(eptrials(time_idx & vid_idx, 9))));
        if hd(wi)<0; hd(wi) = hd(wi)+360; end
    end
    
    if isnan(sum(sum((eptrials(time_idx & vid_idx, 10:17)))))
        obj_vars(wi, :) = nan;
    else
        obj_vars(wi, :) = [mean(eptrials(time_idx & vid_idx, 10:13)) rad2deg(circ_mean(deg2rad(eptrials(time_idx & vid_idx, 14:17))))];
        obj_vars(wi, obj_vars(wi, :)<0) = obj_vars(wi, obj_vars(wi, :)<0) + 360;
    end
    
    %load context (leave nan if window spans multiple contexts)
    current_context = mode(eptrials(time_idx,6));%order 5, id 6
    if length(current_context) == 1
        context_id(wi,1) = current_context;
    end
    
end

end