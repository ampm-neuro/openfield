function [cell_rates, cell_rates_first, cell_rates_second] = ctxt_rates(eptrials, clusters, contexts)
%function ctxt_rates(eptrials, contexts)
%
%calculates the average firing rate for each cell in clusters in each of
%the context conditions in contexts.
%
%also calculates the firing rate for first and second halves in each
%context

%inputs
if size(clusters,1) > size(clusters,2)
    clusters = clusters';
end
if size(contexts,1) > size(contexts,2)
    contexts = contexts';
end

%preallocate
cell_rates = nan(length(clusters), length(contexts));
cell_rates_first = nan(size(cell_rates));
cell_rates_second = nan(size(cell_rates));

%heatmaps_trial(eptrials, clusters(1,1), contexts);

%counters
clust_num = 0;
context_num = 0;

for clust = clusters
    
    %count clusters
    clust_num = clust_num + 1;
    

    for trial_type = contexts
        
        %start time
        starttime = min(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type, 1));
        %end time
        endtime = max(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type, 1));
        
        %timecheck = (endtime - starttime)/60
        %halfway point
        hp = mean([starttime endtime]);
        %hp = min(eptrials(eptrials(:,4)==1 & eptrials(:,5)==trial_type, 1)) + 150;
            
        
        %count contexts
        context_num = context_num + 1;

        %ALL
            %spikes
            %spikes = length(eptrials(eptrials(:,4)==clust & eptrials(:,5)==trial_type & eptrials(:,1)>starttime & eptrials(:,1)<endtime, 1));
            spikes = length(eptrials(eptrials(:,4)==clust & eptrials(:,6)==trial_type & eptrials(:,1)>starttime+30 & eptrials(:,1)<endtime-30, 1));

            %time
            %time = length(eptrials(eptrials(:,4)==1 & eptrials(:,5)==trial_type & eptrials(:,1)>starttime & eptrials(:,1)<endtime, 1))/100;
            time = length(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type & eptrials(:,1)>starttime+30 & eptrials(:,1)<endtime-30, 1))/100;

            %firing rate
            firing_rate = spikes/time;

            %load
            cell_rates(clust_num, context_num) = firing_rate;
            
        %FIRST
            %spikes
            spikes_1 = length(eptrials(eptrials(:,4)==clust & eptrials(:,6)==trial_type & eptrials(:,1)>starttime+30 & eptrials(:,1) < hp, 1));

            %time
            time_1 = length(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type & eptrials(:,1)>starttime+30 & eptrials(:,1) < hp, 1))/100;

            %firing rate
            firing_rate_1 = spikes_1/time_1;

            %load
            cell_rates_first(clust_num, context_num) = firing_rate_1;
            
            
        %SECOND
            %spikes
            spikes_2 = length(eptrials(eptrials(:,4)==clust & eptrials(:,6)==trial_type & eptrials(:,1) > hp & eptrials(:,1)<endtime-30, 1));

            %time
            time_2 = length(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type & eptrials(:,1) > hp & eptrials(:,1)<endtime-30, 1))/100;

            %firing rate
            firing_rate_2 = spikes_2/time_2;

            %load
            cell_rates_second(clust_num, context_num) = firing_rate_2;
    end
    
    %reset context counter
    context_num = 0;
    
end