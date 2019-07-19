function [within, between] = ctxt_rate_diff(eptrials, clusters, context_ids)
%finds the rate difference within and between contexts

%find rates for every cell in clusters
cell_rates = ctxt_rates(eptrials, clusters, context_ids);

%z-score rates
cell_rates = cell_rates - repmat(mean(cell_rates,2), 1,4);
cell_rates = cell_rates ./ repmat(std(cell_rates,[],2), 1,4);

%if its a BW session
if isequal([1; 2; 3; 4], sort(context_ids))
    
   %calculate within and between context difference
   [within, between] = within_between(cell_rates, context_ids, [1 2; 3 4], [1 3; 1 4; 2 3; 2 4]);
    
   %average for output
   within = mean(within)';
   between = mean(between)';
    
%if its IA
elseif ismember([5; 6; 7; 8], sort(context_ids))
    
   %calculate IDENTITY (here, "within") and ARRANGMENT (here, "between") context difference
   [within, between] = within_between(cell_rates, context_ids, [5 6], [7 8]);
   
   %rotate for output
   within = within';
   between = between';
   
   
else
    
   %report error
   error('context session unknown')
   
end

%INTERNAL FUNCTIONS
%
    %find differences
    function rd = rate_diff(cell_rates, context_ids, ctx1, ctx2)
        rd = abs(diff(cell_rates(:,ismember(context_ids, [ctx1 ctx2])), [], 2));        
    end

    %within_between differences
    function [withins, betweens] = within_between(cell_rates, context_ids, within_comparisons, between_comparisons)
       
        %preallocate
        withins = nan(size(within_comparisons,1), size(cell_rates,1));
        betweens = nan(size(between_comparisons,1), size(cell_rates,1));
        
        %iterate through within comparisons
        for wi = 1:size(within_comparisons,1)
            withins(wi, :) = rate_diff(cell_rates, context_ids, within_comparisons(wi,1), within_comparisons(wi,2));
        end
        
        %iterate through between comparisons
        for bi = 1:size(between_comparisons,1)
            betweens(bi, :) = rate_diff(cell_rates, context_ids, between_comparisons(bi,1), between_comparisons(bi,2));
        end
        
    end

end