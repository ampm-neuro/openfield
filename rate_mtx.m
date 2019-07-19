function [rate_matrix, spike_count, spatial_occupancy, rate_matrix_smoothed, HD_proportions] = rate_mtx(eptrials, cluster, bins)
%function [rate_matrix, spike_count, spatial_occupancy] = rate_mtx(eptrials, cluster, bins)
%
% rate_mtx is a simple heatmap plot with adaptive smoothing. It computes
% the firing rate of the neuron 'cluster' within a grid with 'bins' rows and
% 'bins' columns. It depends on precise formatting of the input eptrials,
% where column 1 is time, column 2 is x position, column 3 is y position,
% and column four contains cluster events labeled by their cluster ID, as
% well as timestamps labeled as 1.
%
% rate_mtx outputs the raw rate matrix 'rate_matrix', which is computed by
% divided the output 'spike_count' by the output 'spatial_occupancy'.
%
% rate_mtx also plots a heatmap of the smoothed rate matrix. It uses the
% adaptive smoothing technique first employed by Skagg et al.
%
% ampm 2017

rate_matrix_smoothed = [];
%check inputs
if ~ismember(cluster, unique(eptrials(:,4)))
    input_cluster = cluster
    clusters_in_eptrials = unique(eptrials(:,4))
    warning('Cluster does not appear in eptrials.')
    
    rate_matrix = nan(bins, bins) ;
    spike_count= nan(bins, bins);
    spatial_occupancy = nan(bins, bins);
    HD_proportions = nan(bins, bins, 360);
    rate_matrix_smoothed = nan(bins, bins);
    
    return
    
end

%time minimum
time_min = .5; %s 


%round HDs
eptrials(:,9) = round(eptrials(:,9));
eptrials(eptrials(:,9)==0,9) = 360;

%all spike events in two vectors
xs = eptrials(eptrials(:,4)==cluster, 2);
ys = eptrials(eptrials(:,4)==cluster, 3);

%all time samples in two vectors 
xt = eptrials(eptrials(:,4)==1, 2);
yt = eptrials(eptrials(:,4)==1, 3);



%evenly spaced bins of x and y coordinate ranges
pos_edges_x = linspace(min(eptrials(:,2)), max(eptrials(:,2))+0.000000001, bins+1);
pos_edges_y = linspace(min(eptrials(:,3)), max(eptrials(:,3))+0.000000001, bins+1);

%2d histogram of event counts
spike_count(:,:) = histcounts2(ys, xs, pos_edges_y, pos_edges_x);
spatial_occupancy(:,:) = histcounts2(yt, xt, pos_edges_y, pos_edges_x)./100;

%HD hist counts
for ihd = 1:360
    %all HD samples in two vectors
    xh = eptrials(eptrials(:,9)==ihd, 2);
    yh = eptrials(eptrials(:,9)==ihd, 3);
    HD_proportions(:,:,ihd) = histcounts2(yh, xh, pos_edges_y, pos_edges_x);
end


%flip
spike_count = flipud(spike_count);
spatial_occupancy = flipud(spatial_occupancy);
spatial_occupancy(spatial_occupancy<time_min) = nan;


%divide spikes by time for rate
rate_matrix = spike_count./spatial_occupancy;

%smooth for plotting
%rate_matrix_smoothed = skagg_smooth(spike_count, spatial_occupancy);

%plot
%figure; imagesc(rate_matrix_smoothed); colorbar; title(['Cluster ' num2str(cluster)])


%INTERNAL FUNCTION
%
%adaptive smoothing from Skaggs, McNaughton, Wilson, & Barnes, 1996
    function mtx_smth = skagg_smooth(spike_counts, occupancy_counts)
    % An adaptive smoothing method to optimize the trade-off between blurring 
    % error and sampling error. The firing rate at each bin was estimated
    % by expanding a circle around the point until the radius of the circle 
    % (in bins) is greater than a constant (try 10000) divided by n*sqrt(s),
    % where n is the number of occupancy samples within the circle and s is 
    % the total number of spikes within the circle. With a position sampling 
    % rate of 100 Hz, the firing rate at that point was then set to 100*n*s.
    %
    %Skaggs, McNaughton, Wilson, Barnes 1996, see also Ito 2015. This method is
    %favored by the Mosers.

    %skagg constant
    skg_c = 10000;

    occupancy_counts = occupancy_counts.*100;
    occupancy_counts(occupancy_counts==0) = nan;

    %preallocate rate matrix
    mtx_smth = nan(size(spike_counts));

        %iterate through rows
        for ir = 1:size(spike_counts,1)

            %iterate through columns
            for ic = 1:size(spike_counts,2)

                %skip (leave nan) if no occupancy
                if isnan(occupancy_counts(ir, ic))
                    continue
                end

                %preset test values counter
                radius = 0;
                skagg_val = 1; %arbitrarily higher than radius

                %keep trying until pass skagg test
                while radius < skagg_val

                    %expand smooth mask size
                    radius = radius+1;

                    %Circle with radius centered at ir,ic
                    [mgx, mgy] = meshgrid(1:size(spike_counts,1), 1:size(spike_counts,2));
                    circle_idx = sqrt((mgx-ic).^2+(mgy-ir).^2)<=radius;

                    %sum within circle area
                    circ_spikes = nansum(spike_counts(circle_idx));
                    circ_occupancy = nansum(occupancy_counts(circle_idx));

                    %calculate skagg test
                    skagg_val = skg_c / (circ_occupancy * sqrt(circ_spikes));

                end

                %set smoothed rate
                mtx_smth(ir, ic) = 100 * circ_spikes / circ_occupancy;
            end

        end
    end


end