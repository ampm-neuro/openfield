function [matrices, nanpix] = heatmaps_trial(eptrials, clusts, varargin)
%makes dotplots for all cells clusts
%heatmaps_trial(eptrials, clusts, by trial type?)

%context name index
context_names = {'Black1' 'Black2' 'White1' 'White2' 'Arngmt1' 'Arngmt2' 'Object1' 'Object2'};

matrices = [];

%trial input?
by_trial_type = 0; %no
if nargin == 3
    by_trial_type = 1; %yes 
    trials = varargin{1};
elseif nargin > 3
    error('too many inputs')
end
        
%orient clusts
if size(clusts,1)>size(clusts,2)
    clusts = clusts';
end

%choose plotting method from input
if by_trial_type == 1
    
    %find trial types for plotting
    trial_types = unique(eptrials(~isnan(eptrials(:,6)),6));
    if all(ismember(trials, trial_types))
        trial_types = trials;
    end
    
    %orient trial_types
    if size(trial_types,1)>size(trial_types,2)
        trial_types = trial_types';
    end
    
    %make subplot figure for each cell (Y) by trial types (X)    
    figure
    %h = tight_subplot(3, length(trial_types), [.01 .03],[.1 .01],[.01 .01]);
    if length(clusts)<3
        h = tight_subplot(length(clusts), length(trial_types), [0 0],[0 0],[.01 -.025]);
    else
        h = tight_subplot(3, length(trial_types), [0 0],[0 0],[.01 -.025]);
    end
 
    
    %counters for subplotting
    plotcount = 0;
    count = 0;
    
    %pixles
    pix = 40;
    
    %nanpixles
    nanpix = nan(length(clusts), 4);
    
    %for all cells, all trials
    for clust = clusts
        
        overall_fr = size(eptrials(eptrials(:,4)==clust & ismember(eptrials(:,6), trial_types), 1),1)/(length(eptrials(eptrials(:,4)==1 & ismember(eptrials(:,6), trial_types), 1))/100);
        
        for trial_type = trial_types
            
            %subplot
            count = count+1;
            %tight_subplot(3, length(trial_types), count)
            axes(h(count));
            [~, sc, so, np] = trlfree_heatmap(eptrials(eptrials(:,6)==trial_type, :), clust, pix);
            
            %count nan pixles
            nanpix(clusts == clust, trial_types==trial_type) = np;
            
                matrix_nan_hold = isnan(sc./so) & isinf(sc./so);
                matrix = skagg_smooth(sc, so);
                matrix = smooth2a(matrix,2);
                matrix = inpaint_nans(matrix);
                matrix(matrix_nan_hold) = nan;
            
                matrices = cat(3, matrices, matrix);
                
            imagesc(matrix)
            %axis equal
            %axis tight
            
            %firing rate
            mean_firing_rate = size(eptrials(eptrials(:,4)==clust & eptrials(:,6)==trial_type, 1),1)/(max(eptrials(eptrials(:,6)==trial_type, 1))-min(eptrials(eptrials(:,6)==trial_type, 1)));%(length(eptrials(eptrials(:,4)==1 & eptrials(:,5)==trial_type, 1))/100);

            %titles for subplots
            title([sprintf('%.2f',clust) ' ' context_names{trial_type} ' ' sprintf('%.2f',mean_firing_rate) 'hz'])
            set(gca,'ydir','reverse') %TEMPORARY
            axis off; axis square
            colormap jet
            colorbar
            caxis([0 2*overall_fr])
            %caxis([0 11])
            
            
            %dealing with subplot grouping madness
            if count == length(trials)*3
                plotcount = plotcount + 1;
                print(strcat('dotplot_', num2str(plotcount)),'-depsc')
                count = 0;
                if clust ~= max(clusts)
                    figure
                    %h = tight_subplot(3, length(trial_types), [0 0],[0 0],[.01 -.025]);
                    h = tight_subplot(3, length(trial_types), [0 0],.02,[.01 -.025]);
                end
            end            
        end
    end
    
    
%oh, you didnt want to plot by trial type?    
else

    count = 0;
    for clust = clusts
        count = count+1;
        figure 
        plot(eptrials(:,2), eptrials(:,3), 'Color', [0.8 0.8 0.8])
        hold on
        plot(eptrials(eptrials(:,4)==clust, 2), eptrials(eptrials(:,4)==clust,3), '.', 'Color', 'r')
        axis equal
        title(num2str(clust))
        set(gca,'xdir','reverse')
        %print(strcat('dotplot_', num2str(count)),'-depsc')

    end
    
end