function [vectors_first, vectors_second, vectors_all] = pop_pixlecorr(eptrials, contexts, clusters, bins)
%finds heatmaps for the first ,second halves and full visits to each
%context in contexts for each cell in clusters. 
%
%replaces nans with the mean of the surrounding pixles
%
%indexes (:) to output vectors

%orientation
contexts = (contexts(:))';
clusters = (clusters(:))';

%preallocate 3d matrices (pixles x cells x contexts)
vectors_first = nan(bins^2, length(clusters), length(contexts));
vectors_second = nan(bins^2, length(clusters), length(contexts));
vectors_all = nan(bins^2, length(clusters), length(contexts));

%cell_counter
%celli_count = 0;

%for all cells
for celli = clusters

    %update context counter
    celli_count = find(clusters==celli);
    
    %catch false clusters
    if ~ismember(celli, unique(eptrials(~isnan(eptrials(:,4)),4)))
        
        vectors_first(:, celli_count, :) = zeros(size(vectors_first(:, 1, :)));
        vectors_second(:, celli_count, :) = zeros(size(vectors_first(:, 1, :)));
        vectors_all(:, celli_count, :) = zeros(size(vectors_first(:, 1, :)));
    
        continue
    end
    
    
    %context counter
    contexti_count = 0;
    
    %for all contexts
    for contexti = contexts

        %update context counter
        contexti_count = contexti_count + 1;

        %times
        start_point = min(eptrials(eptrials(:,6)==contexti,1));
        end_point = max(eptrials(eptrials(:,6)==contexti,1));
        mid_point = start_point + (end_point - start_point)/2;

        %matrices
        %{
        matrix_first = hist2_01(eptrials(eptrials(:,1)>=start_point & eptrials(:,1)<mid_point, :), celli, bins);
        matrix_second = hist2_01(eptrials(eptrials(:,1)>=mid_point & eptrials(:,1)<=end_point, :), celli, bins);
        matrix_all = hist2_01(eptrials(eptrials(:,1)>=start_point & eptrials(:,1)<=end_point, :), celli, bins);
        %}
        
        [matrix_first, matrix_first_spike_count, matrix_first_spatial_occupancy]...
            = trlfree_heatmap(eptrials(eptrials(:,1)>=start_point & eptrials(:,1)<mid_point, :), celli, bins, 0);
        matrix_first = smooth_mtx(matrix_first, matrix_first_spike_count, matrix_first_spatial_occupancy);
        
        if isempty(matrix_first)
            break
        end
        
        [matrix_second, matrix_second_spike_count, matrix_second_spatial_occupancy]...
            = trlfree_heatmap(eptrials(eptrials(:,1)>=mid_point & eptrials(:,1)<=end_point, :), celli, bins, 0);
        matrix_second = smooth_mtx(matrix_second, matrix_second_spike_count, matrix_second_spatial_occupancy);
        
        if isempty(matrix_second)
            break
        end
        
        [matrix_all, matrix_all_spike_count, matrix_all_spatial_occupancy]...
            = trlfree_heatmap(eptrials(eptrials(:,1)>=start_point & eptrials(:,1)<=end_point, :), celli, bins, 0);
        matrix_all = smooth_mtx(matrix_all, matrix_all_spike_count, matrix_all_spatial_occupancy);
        
        if isempty(matrix_all)
            break
        end
        

        %fill output matrices
        vectors_first(:, celli_count, contexti_count) = matrix_first(:);
        vectors_second(:, celli_count, contexti_count) = matrix_second(:);
        vectors_all(:, celli_count, contexti_count) = matrix_all(:);

    end
end

%smooth function
    function mtx_out = smooth_mtx(mtx_in, mtx_in_sc, mtx_in_so)
        
        mtx_out = skagg_smooth(mtx_in_sc, mtx_in_so);
        mtx_out = smooth2a(mtx_out, 1);
        mtx_out = inpaint_nans(mtx_out);
        
        
    end



end