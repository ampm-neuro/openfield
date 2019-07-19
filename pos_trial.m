function pos_trial(eptrials, trials)
%makes dotplots for all cells clusts
%heatmaps_trial(eptrials, clusts, by trial type?)

%context name index
context_names = {'Black1' 'Black2' 'White1' 'White2' 'Arngmt1' 'Arngmt2' 'Object1' 'Object2'};

%find trial types for plotting
trial_types = unique(eptrials(~isnan(eptrials(:,6)),6));
if all(ismember(trials, trial_types));
    trial_types = trials;
end

%orient trial_types
if size(trial_types,1)>size(trial_types,2)
    trial_types = trial_types';
end

%make subplot figure for each cell (Y) by trial types (X)    
figure
h = tight_subplot(1, length(trial_types), [0 0],[0 0],[.01 -.025]);

%counters for subplotting
plotcount = 0;
count = 0;

%for all cells, all trials
for trial_type = trial_types

    
    count = count+1;
    subplot(1,length(trials), count)
    
    plot(eptrials(eptrials(:,6)==trial_type,2), eptrials(eptrials(:,6)==trial_type,3), 'Color', [0.8 0.8 0.8])
    axis square

    %titles for subplots
    title(context_names{trial_type})
    set(gca,'ydir','reverse') %TEMPORARY
    %caxis([0 11]) 
    
end
end