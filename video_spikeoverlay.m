function video_spikeoverlay(video_file, eptrials, cluster, context_visit)
%plots video of spikes overlaying video of a context visit

%entrace to context
start_time = min(eptrials(eptrials(:,5)==context_visit,1));

%read video
vid = VideoReader(video_file);

%set video start time
vid.CurrentTime = start_time;


%HOW TO SET BOUNDS? OVERLAY? SAVE?

end

