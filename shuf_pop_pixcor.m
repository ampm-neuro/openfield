function [shuf_means_out_chance, shuf_means_out_ctx, MSEs_all, MSEs_12, MSEs_23] = shuf_pop_pixcor(first_halves, second_halves, shuffs)
%shuffs rates of each cell between pixles and recalculates the correllation
%matrix and corresponding mean correllation

%preallocate
shuf_means_out_chance = nan(shuffs, 3);
shuf_means_out_ctx = nan(shuffs, 3);




%for speed
number_of_contexts = size(first_halves,3);

%iterate shuffs
for si = 1:shuffs

    si
    
    
    %preallocate
    shuf_first_halves_pixles = nan(size(first_halves));
    shuf_first_halves_ctxtids = nan(size(first_halves));
    
    
    %shuf context labels individually for each cell
    for ci = 1:size(first_halves,2)
    
        %shuffle pixles labels
        shuf_first_halves_pixles(:,ci,:) = first_halves(randperm(size(first_halves,1)), ci, :);
        %shuffle context labels
        shuf_first_halves_ctxtids(:,ci,:) = first_halves(:, ci, randperm(number_of_contexts));
    
    end
    
    
    %Compare observed correlation mean to chance (zero)
    %
    % SHOULD SHUFFLE INDEPENDENTLY FOR EACH CELL
    %
    [~, shuf_means_out_chance(si,1), shuf_means_out_chance(si,2), shuf_means_out_chance(si,3)] = all_pop_pixlecorr_plots(shuf_first_halves_pixles, second_halves);
    
    
    %Contrasts
    %
    [~, shuf_means_out_ctx(si,1), shuf_means_out_ctx(si,2), shuf_means_out_ctx(si,3)] = all_pop_pixlecorr_plots(shuf_first_halves_ctxtids, second_halves);

end

%mean squared error
MSEs_all = sum((shuf_means_out_ctx - mean(shuf_means_out_ctx,2)).^2, 2);
MSEs_12 = sum((shuf_means_out_ctx(:,[1 2]) - mean(shuf_means_out_ctx(:,[1 2]),2)).^2, 2);
MSEs_23 = sum((shuf_means_out_ctx(:,[2 3]) - mean(shuf_means_out_ctx(:,[2 3]),2)).^2, 2);













end