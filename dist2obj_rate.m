function [rate_distributions, spike_count_distr, dist_counts, rate_distr_usm] = dist2obj_rate(varargin)
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
%   time_window
%       angle_window is the number of position samples that will be
%       averaged to produce an instantaneous dist from obj. angle_window
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
        time_window = varargin{3};
        smooth_window = 0;
        smooth_multiplier = 0;
    case 4
        eptrials = varargin{1};
        clusters = varargin{2};
        time_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = 1;
    case 5
        eptrials = varargin{1};
        clusters = varargin{2};
        time_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
    case 6
        eptrials = varargin{1};
        clusters = varargin{2};
        time_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
        object_number = varargin{6};
    case 7
        eptrials = varargin{1};
        clusters = varargin{2};
        time_window = varargin{3};
        smooth_window = varargin{4};
        smooth_multiplier = varargin{5};
        object_number = varargin{6};
        ba = varargin{7};
    otherwise
        error(message('Too many input arguments'))
end
%
%force angle_window and smooth_window to be odd
if rem(time_window, 2) == 0
    time_window = ceil(time_window) -1;
end
if rem(smooth_window, 2) == 0
    smooth_window = ceil(smooth_window) -1;
end


w_size = .01;
w_max = .71;

dist_wind_lo = 0:w_size:w_max-w_size;
dist_wind_hi = w_size:w_size:w_max;
dist_wind_mid = mean([dist_wind_lo(1); dist_wind_hi(1)]):...
    w_size:mean([dist_wind_lo(end); dist_wind_hi(end)]);
%ba = [0 1 0 1];


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
window_length = time_window * window_to_seconds;
%
%cluster vectors
clust_pop = unique(eptrials(eptrials(:,4)>1,4));
desired_clusters = sort(clusters);
%
%time point vectors
all_time_pts = eptrials(eptrials(:,4)==1, 1);
mid_time_pts = all_time_pts((floor(time_window/2)+1):(length(all_time_pts) - floor(time_window/2)));
%
%calculate instantaneous dist from obj
inst_dist = eptrials(eptrials(:,4) == 1, object_number);
inst_dist((length(inst_dist) - floor(time_window/2)+1):end) = [];%trim end
inst_dist(1:floor(time_window/2)) = [];%trim begining
inst_dist = inst_dist(1:iteration_jump:end);

%%%TEST TEST TEST
%{
mid_time_pts(inst_dist<dist_wind_lo(1) | inst_dist>dist_wind_hi(end)) = [];
inst_dist(inst_dist<dist_wind_lo(1) | inst_dist>dist_wind_hi(end)) = [];
%}

for i = 1:length(dist_wind_lo)
    %renumber inst_hd to rounded values of hd_angles
    inst_dist(inst_dist>=dist_wind_lo(i) & inst_dist<=dist_wind_hi(i)) = dist_wind_mid(i);
end

%
%preallocate spike counts vector
spike_counts = nan(length(clust_pop), length(mid_time_pts));%length(all_time_pts)
%
%preallocate rate distribution matrix
spike_count_distr = nan(length(clust_pop), length(dist_wind_mid));
rate_distr = nan(length(clust_pop), length(dist_wind_mid));
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
    wndw_index = eptrials(:,1) >= all_time_pts(time_sample)...
        & eptrials(:,1) <= all_time_pts(time_sample + time_window -1);
    
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
inst_dist = inst_dist(mode(spike_counts~=.5,1));
spike_counts = spike_counts(:, mode(spike_counts~=.5,1));

%USE SPIKE COUNTS TO BUILD RATE DISTRIBUTION MATRIX
%
%number of time windows corresponding to each head direction 1:360
dist_counts = histc(inst_dist, dist_wind_mid);
if sum(dist_counts==0)>1
    warning('missing hd values; data ommitted')
    dist_counts(dist_counts==0) = nan;
end

%number of spikes corresponding to each distance (bin) from obj
for id = 1:length(dist_wind_mid)
    dist_from_obj = dist_wind_mid(id);
    
    %avoiding empty sets
    if sum(inst_dist == dist_from_obj) > 0
        %rate distribution at degree dir is equal to the sum of all spike
        %counts occuring in windows corresponding to that head direction.
        spike_count_distr(:, id) = sum(spike_counts(:, inst_dist == dist_from_obj), 2);
    else
        spike_count_distr(:, id) = zeros(size(rate_distr(:, id)));
    end
end
%
%time spent at each head direction
dist_time = dist_counts.*window_length;

%divide counts by time to get rates
rate_distr = bsxfun(@rdivide, spike_count_distr, dist_time');
rate_distr_usm = rate_distr;

%CIRCULAR SMOOTH RATE DISTRIBUTION MATRIX
%
%if inputs call for smoothing
if smooth_multiplier > 0
    
    %smooth one row (cluster) at a time
    for current_cluster = 1:size(rate_distr,1)
    
        %smooth_multiplier is the number of times smoothing is performed
        for smoothing_round = 1:smooth_multiplier
            rate_distr(current_cluster,:) = nanfastsmooth(rate_distr(current_cluster,:), smooth_window);
        end 
    end
end

%isolate distributions from desired cells
rate_distributions = rate_distr(ismember(clust_pop, desired_clusters),:);

end