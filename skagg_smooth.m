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
            
            if radius > size(spike_counts,2)/2
                mtx_smth = [];
                return
            end
                

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