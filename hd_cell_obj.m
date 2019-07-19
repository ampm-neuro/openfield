function [rate_distributions, spike_count_distr, dist_counts, rate_distr_usm] = hd_cell_obj(varargin)
% function hd_cell(eptrials, clusters, angle_window, smooth_window, smooth_multiplier)
%
% hd_cell plots the firing rate of a single cells "clusters" over 360 
% degrees of head direction. eptrials, clusters, and angle_window are 
% required inputs, while smooth_window and smooth_multiplier are optional.
%
% INPUTS:
%   eptrials
%       eptrials is a data preprocessing matrix output by the function
%       "trials"
%   clusters
%       clusters is a vector of ID numbers of cells that will be plotted
%   angle_window
%       angle_window is the number of head direction samples that will be
%       averaged to produce an instantaneous head direction. angle_window
%       must be odd, and will be reduced by one if an even input is given. 
%       Try: 11.
%   smooth_window
%       smooth_window is an optional input that defines the width of the
%       sliding average used to smooth the distribution of firing rates
%       over 360 degrees. If only 3 inputs, no smoothing will occur.
%   smooth_multiplier
%       smooth_multiplier is an optional input that defines the number of
%       times smoothing will be repeated. If only 4 inputs, 
%       smooth_multiplier defaults to 1.
%        
w_size = 5;
hd_angles_lo = 0:w_size:360-w_size;
hd_angles_hi = w_size:w_size:360;
hd_angles_mid = mean([hd_angles_lo(1); hd_angles_hi(1)]):w_size:mean([hd_angles_lo(end); hd_angles_hi(end)]);
ba = [0 1];
%CHECK INPUTS
%
switch nargin
    case 0
        error(message('Need more input arguments'))
    case 1
        error(message('Need more input arguments'))
    case 2
        error(message('Need more input arguments'))
    case 3
        eptrials = varargin{1};
        clusters = varargin{2};
        angle_window = varargin{3};
        smooth_window = 0;
        smooth_multiplier = 0;
    case 4
        eptrials = varargin{1};
        clusters = varargin{2};
        angle_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = 1;
    case 5
        eptrials = varargin{1};
        clusters = varargin{2};
        angle_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
    case 6
        eptrials = varargin{1};
        clusters = varargin{2};
        angle_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
        object_number = varargin{6};
    case 7
        eptrials = varargin{1};
        clusters = varargin{2};
        angle_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
        object_number = varargin{6};
        ba = varargin{7};
    otherwise
        error(message('Too many input arguments'))
end
%
%force angle_window and smooth_window to be odd
if rem(angle_window, 2) == 0
    angle_window = ceil(angle_window) -1;
end
if rem(smooth_window, 2) == 0
    smooth_window = ceil(smooth_window) -1;
end


%bound arena
%can include only subsets of positions (e.g. remove walls and corners)
%ba = [0 1]; %low hi

%WORKING VARIABLES
%
%how many video samples to skip in each iteration
iteration_jump = 5;
%
%window length
window_to_seconds = 0.01; %can change time conversion here
window_length = angle_window * window_to_seconds;
%
%cluster vectors
clust_pop = unique(eptrials(eptrials(:,4)>1,4));
desired_clusters = sort(clusters);
%
%time point vectors
all_time_pts = eptrials(eptrials(:,4)==1, 1);
%begin_time_pts = all_time_pts(1:floor(angle_window/2));
%end_time_pts = all_time_pts((length(all_time_pts) - floor(angle_window/2)+1):end);
mid_time_pts = all_time_pts((floor(angle_window/2)+1):(length(all_time_pts) - floor(angle_window/2)));
%
%calculate instantaneous head direction
inst_hd = circ_smooth(eptrials(eptrials(:,4) == 1, object_number), [0 360], angle_window);
inst_hd((length(inst_hd) - floor(angle_window/2)+1):end) = [];%trim end
inst_hd(1:floor(angle_window/2)) = [];%trim begining
inst_hd = inst_hd(1:iteration_jump:end);
for i = 1:length(hd_angles_lo)
    inst_hd(inst_hd>=hd_angles_lo(i) & inst_hd<=hd_angles_hi(i)) = hd_angles_mid(i);
end

%renumber inst_hd to rounded values of hd_angles

%
%preallocate spike counts vector
spike_counts = nan(length(clust_pop), length(mid_time_pts));%length(all_time_pts)
%
%preallocate rate distribution matrix
spike_count_distr = nan(length(clust_pop), length(hd_angles_mid));
rate_distr = nan(length(clust_pop), length(hd_angles_mid));
%
%set counter for angle window loop
iteration = 0;%length(begin_time_pts);


%CALCULATE INSTANTANEOUS SPIKE COUNTS
%
%iterate through middle time points (ignoring very ends)
for time_sample = 1:iteration_jump:length(mid_time_pts)
    
    %counter
    iteration = iteration + 1;

    %current time window
    %wndw_low = eptrials(eptrials(:,4) == 1 & eptrials(:,1) == all_time_pts(time_sample), 1);
    %wndw_high = eptrials(eptrials(:,4) == 1 & eptrials(:,1) == all_time_pts(time_sample + angle_window -1), 1);
    wndw_index = eptrials(:,1) >= all_time_pts(time_sample)...
        & eptrials(:,1) <= all_time_pts(time_sample + angle_window -1);
    
    if mean(eptrials(wndw_index,2))>ba(1) && mean(eptrials(wndw_index,2))<ba(2) && mean(eptrials(wndw_index,3))>ba(3) && mean(eptrials(wndw_index,3))<ba(4)
        %identify spike counts during window for each cluster
        current_spk_evts = eptrials(wndw_index, 4);
        spike_counts(:,iteration) = histc(current_spk_evts(~isnan(current_spk_evts)), unique(clust_pop)); 
    else
        spike_counts(:,iteration) = repmat(.5, size(spike_counts(:,1))); %arbitrary impossible spike count to flag bin
    end
end
spike_counts = spike_counts(:, ~isnan(sum(spike_counts, 1)));

%down sample inst_hd to only match timewindows where I have spike counts
%this removes excluded windows
inst_hd = inst_hd(mode(spike_counts~=.5,1));
spike_counts = spike_counts(:, mode(spike_counts~=.5,1));

%USE SPIKE COUNTS TO BUILD RATE DISTRIBUTION MATRIX
%
%number of time windows corresponding to each head direction 1:360
dist_counts = histc(inst_hd, hd_angles_mid);
if sum(dist_counts==0)>1
    warning('missing hd values; data ommitted')
    dist_counts(dist_counts==0) = nan;
end

%number of spikes corresponding to each head direction 1:360
for id = 1:length(hd_angles_mid)
    dir = hd_angles_mid(id);
    
    %avoiding empty sets
    if sum(double(inst_hd == dir),2) > 0
    
        %rate distribution at degree dir is equal to the sum of all spike
        %counts occuring in windows corresponding to that head direction.
        spike_count_distr(:, id) = sum(spike_counts(:, inst_hd == dir), 2);
    else
        spike_count_distr(:, id) = zeros(size(rate_distr(:, id)));
    end
end
%
%time spent at each head direction
hd_time = dist_counts.*window_length;

%divide counts by time to get rates
rate_distr = bsxfun(@rdivide, spike_count_distr, hd_time);
rate_distr_usm = rate_distr;

%CIRCULAR SMOOTH RATE DISTRIBUTION MATRIX
%
%if inputs call for smoothing
if smooth_multiplier > 0
    %pad rate distribution matrix
    rate_distr = [rate_distr rate_distr rate_distr];
    
    %
    %smooth one row (cluster) at a time
    for current_cluster = 1:size(rate_distr,1)
    
        %smooth_multiplier is the number of times smoothing is performed
        for smoothing_round = 1:smooth_multiplier
            %rate_distr(current_cluster,:) = circ_smooth(rate_distr(current_cluster,:), [0 360], smooth_window );
            %rate_distr(current_cluster,:) = smooth(rate_distr(current_cluster,:), smooth_window);
            rate_distr(current_cluster,:) = nanfastsmooth(rate_distr(current_cluster,:), smooth_window);
        end 
    end
    %
    %extract the now-smoothed matrix from the padded matrix
    rate_distr = rate_distr(:,((size(rate_distr,2)/3)+1):(size(rate_distr,2)*(2/3)));
    
end
%}
%clust_pop = clust_pop
%desired_clusters = desired_clusters

%isolate distributions from desired cells
rate_distributions = rate_distr(ismember(clust_pop, desired_clusters),:);

%PLOT DESIRED RATE DISTRIBUTIONS
%
info_content = nan(length(clusters),1);

num_figs = size(findobj('Type', 'figure'),1);
target_figs = (num_figs-length(clusters)+1):num_figs;
%{
for cluster_number = 1:length(clusters)
    target_fig = target_figs(cluster_number);
    
    %active clusters
    current_cluster = desired_clusters(cluster_number);
    rate_distr_row = find(clust_pop==current_cluster, 1);

    %plot
    %{
    figure(target_fig)
    hold on
    if ismember(mode(eptrials(:,6)), [1 5])
        plot(rate_distr(rate_distr_row,:), 'color', [0    0.4470    0.7410], 'linewidth', 1.1)
    elseif ismember(mode(eptrials(:,6)), [2 6])
        plot(rate_distr(rate_distr_row,:), 'color', [0 .24 .54], 'linewidth', 1.1)    
    elseif ismember(mode(eptrials(:,6)), [3 7])
        plot(rate_distr(rate_distr_row,:), 'color', [0.8500    0.3250    0.0980], 'linewidth', 1.1) 
    elseif ismember(mode(eptrials(:,6)), [4 8])
        plot(rate_distr(rate_distr_row,:), 'color', [0.6500    0.1250    0], 'linewidth', 1.1) 
    else
        plot(rate_distr(rate_distr_row,:), 'color', [0.7 0.7 0.7], 'linewidth', 1.1)
    end
    set(gca,'TickLength',[0, 0]); box off
    %axis([0 360 0 max(rate_distr(rate_distr_row,:))*1.2])
    set(gca,'XTick', [1 90 180 270 360], 'fontsize', 17)
    title(num2str(current_cluster), 'fontsize', 20)
    ylabel('Mean Firing Rate (Hz)', 'fontsize', 20)
    xlabel('Head Direction (clockwise degrees)', 'fontsize', 20) 
    %}
    
    %probability of facing each direction
    %(occupying each direction bin)
    direction_p = hd_counts./nansum(hd_counts(:));
    
    %remove unvisited bins
    rd_hold = rate_distributions(cluster_number, direction_p>0);
    p_hold = direction_p(direction_p>0);
    
    %overall mean firing rate
    R = nansum(nansum(rd_hold.*p_hold));
    
    %for each direction bin
    info_content_dir = nan(360,1);

    %all directions
    ad = 1:360;
    
    %info content for each bin
    for i_dir = ad(direction_p>0)
        Pi = direction_p(i_dir);%probability of occupancy of bin i
        Ri = rate_distributions(cluster_number, i_dir);%mean firing rate for bin i
        if Ri==0
            Ri = 0.0000001;
        end
        rr = Ri/R;%relative rate
        info_content_dir(i_dir) = Pi*(rr)*log2(rr);%load info content for each bin

    end
    info_content(cluster_number) = sum(info_content_dir(direction_p>0));

        
end
%}



end