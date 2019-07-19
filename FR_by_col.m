function [rates_out, col_out] = FR_by_col(eptrials, clusters, col, bins, tw_size, slide_size, varargin)
%plots firing rate of cell 'cluster' over equal sized bins, of number
%'bins', of items in column 'col'
%
% for head direction cell plots use col 9
%
% varargin
% 2 item vector stating bin range. for HD use [0 360]

%varargin
if nargin > 6
    range = varargin{1};
else
    range =[min(eptrials(:,col)) max(eptrials(:,col))]; 
end

%COUNT SPIKES WITH A SLIDING WINDOW
%
%sliding window
%tw_size = 500; %ms
%slide_size = 500; %ms

%PREALLOCATE
col_val = nan((sum(eptrials(:,4)==1)), 1);
spike_counts = nan((sum(eptrials(:,4)==1)), length(clusters));

%all times
all_times = eptrials(eptrials(:,4)==1,1);

    %preindex
    vid_idx = eptrials(:,4)==1;
for tw = 1:(slide_size/10):length(eptrials(eptrials(:,4)==1,1))-(tw_size/10)
    
    %tw index
    tw_idx = eptrials(:,1)>all_times(tw) & eptrials(:,1)<=all_times(tw+(tw_size/10));
    
    %current col_val
    col_val(tw,1) = mean(eptrials(tw_idx & vid_idx, col));

    %count spikes in window (sets upper bin past all clusters)
    dummy_clusters = [clusters; clusters(end)+.0001];
    spike_counts(tw,:) = histcounts(eptrials(tw_idx & ~vid_idx,4), dummy_clusters);
end

% ~isnan idx
nnan_idx = ~isnan(spike_counts(:,1)) & ~isnan(col_val(:,1));
spike_counts = spike_counts(nnan_idx,1:length(clusters));
col_val = col_val(nnan_idx);

%CALCULATE FIRING RATES AT EACH BIN OF COL_VAL       
%
%bin edges
edges = linspace(range(1), range(2), bins+1);

%count the number of col_vals that fall in each bin
[cv_histc, ~, histcount_idx] = histcounts(col_val, edges);
cv_histc = cv_histc';

%preallocate spike count matrix
sc_histc = nan(length(cv_histc), length(clusters));

%count the number of spikes that fall in each colval bin
for ie = 1:length(cv_histc)
   sc_histc(ie, :) = sum(spike_counts(histcount_idx==ie, :));
end

%divide the spike counts by the time spent in each bin
time_in_each_bin = (repmat(cv_histc,1,length(clusters)).*(tw_size/1000));
    min_time = 2; %seconds
    time_in_each_bin(time_in_each_bin < min_time) = NaN;
    sc_histc = sc_histc./time_in_each_bin;


for ic = 1:size(sc_histc,2)
    rates_out(:,ic) = sc_histc(:,ic);
end
col_out = linspace(mean([edges(1) edges(2)]),  mean([edges(end-1) edges(end)]), size(rates_out,1));
%mean([edges(1) edges(2)]):abs(diff([edges(1) edges(2)])):mean([edges(end-1) edges(end)]);
col_out = col_out';



end