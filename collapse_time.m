function time_collapsed_rates = collapse_time(all_counts, context_idx)
%takes ouput of all_countwindow (all_counts) and collapses across the
%time dimension (the vector through the first and last iti). output can be
%used for distance measures or can be pca'd and then can be plotted. 
%
%for all_counts use "all_counts_bw_3" or similar
%for context_idx use "context_bw_3" or similar

all_counts = all_counts(:, ~isnan(all_counts(1,:)));


%collapse
%
time_collapsed_rates = all_counts*null(mean(all_counts(ismember(context_idx, 12),:)) - mean(all_counts(ismember(context_idx, 10),:)));


%plot
%{
    tcr_pca = pca(time_collapsed_rates');

    %color
    colors = [17 17 193;...
              38 166 230;...
              202 16 16;...
              228 146 40;...
              50 50 50;...
              115 115 115;...
              154 154 154]./255;


    figure; hold on
    count = 0;
        for ctxti = unique(context_idx)'
            count = count+1;
            plot3(tcr_pca(context_idx==ctxti,1), tcr_pca(context_idx==ctxti,2), tcr_pca(context_idx==ctxti,3),'.', 'MarkerSize', 8, 'Color', colors(count,:)); 
            h{count} = plot3(mean(tcr_pca(context_idx==ctxti,1)), mean(tcr_pca(context_idx==ctxti,2)), mean(tcr_pca(context_idx==ctxti,3)),'.', 'MarkerSize', 60, 'Color', colors(count,:)); 
        end
    xlabel('pc1')
    ylabel('pc2')
    zlabel('pc3')
    legend show
    if ismember(1, context_idx)
        legend([h{1}, h{2}, h{3}, h{4}, h{5}, h{6}, h{7}], 'Black1', 'Black2', 'White1', 'White2', 'ITI1', 'ITI2', 'ITI3')
    else
        legend([h{1}, h{2}, h{3}, h{4}, h{5}, h{6}, h{7}], 'Arng1', 'Arng2', 'Obj1', 'Obj2', 'ITI1', 'ITI2', 'ITI3')
    end
    axis square
end
%}