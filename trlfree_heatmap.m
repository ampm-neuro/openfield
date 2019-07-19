function [rate_matrix, spike_count, spatial_occupancy, p_nan] = trlfree_heatmap(eptrials, clust, bins, varargin)
%makes heatmap
%
%eptrials, clust, bins (same for x and y), and varargin figure_on (default
%is off
%
%plot_figure should be 1 for a heatmap, 0 otherwise
%min_max is [minX, maxX; minY, maxY]

%INPUT
%
if nargin == 4
    figure_on = varargin{1};%plot figure
elseif nargin == 5
    figure_on = varargin{1};%plot figure
    object_condition = varargin{2};
elseif nargin == 6
    figure_on = varargin{1};%plot figure
    object_condition = varargin{2};
    min_max = varargin{3}; %input min and max eptrials values
    
else
    figure_on = 0;
end


%COMPUTE OUTPUT
%
%all spike events in two vectors
xs = eptrials(eptrials(:,4)==clust, 2);
ys = eptrials(eptrials(:,4)==clust, 3);

%all time samples in two vectors 
xt = eptrials(eptrials(:,4)==1, 2);
yt = eptrials(eptrials(:,4)==1, 3);

%evenly spaced bins of x and y coordinate ranges
if exist('min_max', 'var')
    pos_edges_x = linspace(min_max(1,1), min_max(1,2)+0.000000001, bins+1);
    pos_edges_y = linspace(min_max(2,1), min_max(2,2)+0.000000001, bins+1);
else
    pos_edges_x = linspace(min(eptrials(:,2)), max(eptrials(:,2))+0.000000001, bins+1);
    pos_edges_y = linspace(min(eptrials(:,3)), max(eptrials(:,3))+0.000000001, bins+1);
end

%2d histogram of event counts
spike_count = histcounts2(ys, xs, pos_edges_y, pos_edges_x);
spatial_occupancy = histcounts2(yt, xt, pos_edges_y, pos_edges_x)./100;

%flip
spike_count = flipud(spike_count);
spatial_occupancy = flipud(spatial_occupancy);

%divide spikes by time for rate
rate_matrix = spike_count./spatial_occupancy;

%proportion of unvisited pixles
p_nan = sum(sum(isnan(rate_matrix)))/ numel(rate_matrix);


%PLOT
%
%plot if desired
if figure_on == 1
    
    %smooth
    rate_matrix_smoothed = skagg_smooth(spike_count, spatial_occupancy);
    rate_matrix_smoothed = inpaint_nans(rate_matrix_smoothed);
    
    %zscore
    rate_matrix_smoothed = rate_matrix_smoothed - ...
        repmat(mean(rate_matrix_smoothed(:)), size(rate_matrix_smoothed));
    rate_matrix_smoothed = rate_matrix_smoothed./std(rate_matrix_smoothed(:));
    
    %figure; 
    imagesc(rot90(rate_matrix_smoothed,2)); %colorbar; title rate
    
    if exist('object_condition', 'var')
        hold on; object_grid(object_condition, bins);
    end
end

end