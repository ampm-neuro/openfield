%
time_window =  .25; 

%all_counts vars from 
%{
[all_counts_bw, all_countsz_ia, context_bw, context_ia,...
all_counts_bw_z, all_countsz_ia_z, region_context_bw, region_context_ia,...
all_pos_bw, all_hd_bw, all_pos_ia, all_visit_order_bw, all_visit_order_ia] = ...
ALL_countwindow(time_window);
%}

%distances
%CONSIDER USING A CLOUD TO CLOUD DISTANCE
%OR DIVIDING BY MEAN WITHIN SD
method = [1, 2];
tws = (12*60) / time_window;

%
mtx = all_counts_bw;
%mtx = mtx(:, region_context_bw == 2); %regional differences 1 = RSA, 2 = RGB
%mtx = mtx(:, ~isnan(mtx(1,:)));

%time collapse
%mtx = collapse_time(mtx, context_bw);
%mtx = collapse_time_individual(mtx, context_bw, all_visit_order_bw);

%mtx = mtx(~ismember(context_bw, [10 11 12]),:);
%context_bw = context_bw(~ismember(context_bw, [10 11 12]),:);

mtx = zscore_mtx(mtx);


%PCA FIG
%{
colors = [17 17 193;38 166 230;202 16 16;228 146 40;50 50 50;115 115 115;154 154 154]./255;
[mtx_pca, ~, ~, ~, explained_bw_pca] = pca(mtx');
%plot BW
figure; hold on
count = 0;
for ctxti = unique(context_bw)'
    count = count+1;
    plot3(mtx_pca(context_bw==ctxti,1), mtx_pca(context_bw==ctxti,2), mtx_pca(context_bw==ctxti,3),'.', 'MarkerSize', 8, 'Color', colors(count,:));
    plot3(mean(mtx_pca(context_bw==ctxti,1)), mean(mtx_pca(context_bw==ctxti,2)), mean(mtx_pca(context_bw==ctxti,3)),'.', 'MarkerSize', 60, 'Color', colors(count,:)); 
end
%}


%[w_bar, b_bar] = spkct_distances_bw(mtx, context_bw, [1 2; 3 4], [1 3; 2 4; 2 3; 1 4], method(1), method(2));


%classifier first visits are used to train, second to test
%
ctx_context_bw = context_bw(ismember(context_bw, 1:4), :);
ctx_context_bw(ctx_context_bw==2) = 1;
ctx_context_bw(ismember(ctx_context_bw, [3 4])) = 2;

data_train_bw_1 = mtx(ismember(context_bw, [1 3]),:);
data_test_bw_1 = mtx(ismember(context_bw, [2 4]),:);
class_train_bw_1 = ctx_context_bw(ismember(context_bw, [1 3]),:);
class_test_bw_1 = ctx_context_bw(ismember(context_bw, [2 4]),:);
[p_correct_bw_1, assignments_bw_1, distances_bw_1, class_test_bw_1] = md_class(data_train_bw_1, class_train_bw_1, data_test_bw_1, class_test_bw_1);
title 1

data_train_bw_2 = mtx(ismember(context_bw, [2 4]),:);
data_test_bw_2 = mtx(ismember(context_bw, [1 3]),:);
class_train_bw_2 = ctx_context_bw(ismember(context_bw, [2 4]),:);
class_test_bw_2 = ctx_context_bw(ismember(context_bw, [1 3]),:);
[p_correct_bw_2, assignments_bw_2, distances_bw_2, class_test_bw_2] = md_class(data_train_bw_2, class_train_bw_2, data_test_bw_2, class_test_bw_2);
title 2
%}

%combined plot
%plot_bw_mdclass(distances_bw_1, distances_bw_2, class_test_1, class_test_2);

%cell_e_bw{1} = [distances_bw(ctx_context_bw==1, 1); distances_bw(ctx_context_bw==2, 2)];
%cell_e_bw{2} = [distances_bw(ctx_context_bw==1, 2); distances_bw(ctx_context_bw==2, 1)];
%errorbar_barplot(cell_e_bw)
%[~,~,~, t_stats_bw] = ttest(cell_e_bw{1}, cell_e_bw{2});
%}
%

%classification errors
%{
b_errors = sum(ctx_context_bw == 1 & assignments_bw == 2) + sum(ctx_context_bw == 2 & assignments_bw == 1);
w_errors = sum(ctx_context_bw == 3 & assignments_bw == 4) + sum(ctx_context_bw == 4 & assignments_bw == 3);
within_errors_bw = b_errors + w_errors;
between_errors_bw = sum(assignments_bw~=ctx_context_bw) - within_errors_bw;
figure; bar([within_errors_bw between_errors_bw]./(within_errors_bw + between_errors_bw)); xticks([1 2]); xticklabels({'Within', 'Between'}); ylabel('Number of errors')
%}
%}
%}
%
mtx = all_countsz_ia;

mtx = collapse_time(mtx, context_ia);
%mtx = collapse_time_individual_ia(mtx, context_ia, all_visit_order_ia);

%mtx = mtx(~ismember(context_ia, [10 11 12]),:);
%context_ia = context_ia(~ismember(context_ia, [10 11 12]),:);



mtx = zscore_mtx(mtx);
%mtx = mtx(:, region_context_bw == 2); %regional differences 1 = RSA, 2 = RGB
%mtx = mtx(:, ~isnan(mtx(1,:)));


%PCA FIG
%
colors = [17 17 193;38 166 230;202 16 16;228 146 40;50 50 50;115 115 115;154 154 154]./255;
[mtx_pca, ~, ~, ~, explained_ia_pca] = pca(mtx');
%plot BW
figure; hold on
count = 0;
for ctxti = unique(context_ia)'
    count = count+1;
    plot3(mtx_pca(context_ia==ctxti,1), mtx_pca(context_ia==ctxti,2), mtx_pca(context_ia==ctxti,3),'.', 'MarkerSize', 8, 'Color', colors(count,:));
    plot3(mean(mtx_pca(context_ia==ctxti,1)), mean(mtx_pca(context_ia==ctxti,2)), mean(mtx_pca(context_ia==ctxti,3)),'.', 'MarkerSize', 60, 'Color', colors(count,:)); 
end
%}


[~, ~, a_bar, o_bar, ao_bar] = spkct_distances_ia(mtx, context_ia, [5; 6; 7; 8], [5 6; 7 8; 5 7; 6 8; 5 8; 6 7], method(1), method(2));

%}

%classifier
data_train_ia_1 = [];
data_test_ia_1 = [];
class_train_ia_1 = [];
class_test_ia_1 = [];

data_train_ia_2 = [];
data_test_ia_2 = [];
class_train_ia_2 = [];
class_test_ia_2 = [];

mtx = mtx_pca(:,1:3);

for ctx_ia = 5:8
    
    all_data_hold = mtx(ismember(context_ia, ctx_ia),:);
    firsthalf_data_hold = all_data_hold(1:floor(size(all_data_hold,1)/2),:);
    secondhalf_data_hold = all_data_hold(ceil(size(all_data_hold,1)/2):end,:);   
    
    data_train_ia_1 = [data_train_ia_1; firsthalf_data_hold];
    data_test_ia_1 = [data_test_ia_1; secondhalf_data_hold];
    class_train_ia_1 = [class_train_ia_1; repmat(ctx_ia, size(firsthalf_data_hold, 1), 1)];
    class_test_ia_1 = [class_test_ia_1; repmat(ctx_ia, size(secondhalf_data_hold, 1), 1)];

    data_train_ia_2 = [data_train_ia_2; secondhalf_data_hold];
    data_test_ia_2 = [data_test_ia_2; firsthalf_data_hold];
    class_train_ia_2 = [class_train_ia_2; repmat(ctx_ia, size(secondhalf_data_hold, 1), 1)];
    class_test_ia_2 = [class_test_ia_2; repmat(ctx_ia, size(firsthalf_data_hold, 1), 1)];
    
end


[p_correct_ia_1, assignments_ia_1, distances_ia_1, class_test_ia_1] = md_class(data_train_ia_1, class_train_ia_1, data_test_ia_1, class_test_ia_1);
title 1

[p_correct_ia_2, assignments_ia_2, distances_ia_2, class_test_ia_2] = md_class(data_train_ia_2, class_train_ia_2, data_test_ia_2, class_test_ia_2);
title 2




%ctx_mtx_ia = mtx(ismember(context_ia, 5:8), :);
%ctx_context_ia = context_ia(ismember(context_ia, 5:8));
%[p_correct_ia, assignments_ia, distances_ia] = md_classify(ctx_mtx_ia, ctx_context_ia);
%plot_ia_mdclass_dotplot

%classification errors ia
a_errors = sum(class_test_ia_1 == 5 & assignments_ia_1 == 6) + sum(class_test_ia_1 == 6 & assignments_ia_1 == 5)...
    + sum(class_test_ia_2 == 5 & assignments_ia_2 == 6) + sum(class_test_ia_2 == 6 & assignments_ia_2 == 5);
i_errors = sum(class_test_ia_1 == 7 & assignments_ia_1 == 8) + sum(class_test_ia_1 == 8 & assignments_ia_1 == 7)...
    + sum(class_test_ia_2 == 7 & assignments_ia_2 == 8) + sum(class_test_ia_2 == 8 & assignments_ia_2 == 7);
%ia_errors = sum(assignments_ia~=ctx_context_ia) - a_errors - i_errors;
%figure; bar([i_errors a_errors ia_errors]./((1 - p_correct_ia)*length(ctx_context_ia))); xticks([1 2]); xticklabels({'Object', 'Arrangement'}); ylabel('Number of errors')
%figure; bar([i_errors a_errors]./((1 - p_correct_ia)*length(ctx_context_ia))); xticks([1 2]); xticklabels({'Object', 'Arrangement'}); ylabel('Number of errors')


error yo

%reshape and average
%b_bar = reshape(b_bar, tws, 4); b_bar = mean(b_bar,2); cell_e{1} = b_bar;
%w_bar = reshape(w_bar, tws, 4);w_bar = mean(w_bar,2); cell_e{2} = w_bar;
a_bar = reshape(a_bar, tws, 2); a_bar = mean(a_bar,2); cell_e{3} = a_bar;
o_bar = reshape(o_bar, tws, 2); o_bar = mean(o_bar,2); cell_e{4} = o_bar;
ao_bar = reshape(ao_bar, tws, 8); ao_bar = mean(ao_bar,2); cell_e{5} = ao_bar;


%figure
%
errorbar_barplot(cell_e)
%ylim([1 1+(mean(b_bar)-1)*1.15])
xlim([.25 5.75])
xticklabels({'Between', 'Within', 'Arrange', 'Object', 'Arg&Obj'})
xlabel('Context Comparison')
ylabel('Distance (Standard Deviations)')
%}

%single cell distribution
[rate_cell_out, space_cell_out] = ALL_singlcell_ctxtdiscrm_distr;
figure; histogram(rate_cell_out{1}, 60, 'normalization', 'probability');
set(gca,'TickLength',[0, 0]); box off
xlim([-30 30])
axis square
xlabel('zscore difference')
ylabel('proportion of cells')
title('BW context discrimination')
%}

function [p_correct_bw, assignments_bw, distances_bw, class_test_bw] = md_class(data_train_bw, class_train_bw, data_test_bw, class_test_bw)
[p_correct_bw, assignments_bw, distances_bw, class_test_bw] = md_classify_traintest(data_train_bw, class_train_bw, data_test_bw, class_test_bw);
%plot_bw_mdclass_dotplot_12_halves
distances_bw = distances_bw;
figure; bar([sum(assignments_bw(class_test_bw==1)==1)/sum(class_test_bw==1) sum(assignments_bw(class_test_bw==1)==2)/sum(class_test_bw==1)...
    sum(assignments_bw(class_test_bw==2)==2)/sum(class_test_bw==2) sum(assignments_bw(class_test_bw==2)==1)/sum(class_test_bw==2)])
end