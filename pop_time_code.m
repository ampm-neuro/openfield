%plots evidence for a rate code for each visit 1-4
%takes the first 6min and the last 6min of every visit (visits should be 12
%min long, but are often slightly longer).
%sessions with insufficient pre or post visit time (below) are excluded

%amount of time for rate windows
time_window =  .250; 

%time to include before first visit. This sets reference time point for 
%populations distances (important that it's far enough from start of v1)
pre_visit1_time = 60;%seconds 


%time to include from ITIs. Half of the time from the start of every iti,
%half from the end of every iti.
iti_time = 180;

%time to include after last visit.
post_visit4_time = 60;

%full drift code
%easy set to just BW or just IA or both, but you have to go into the
%code... dum dum dum.
%
[rm, rm_z, ctxid, visit_idx, half_idx] = ALL_countwindow_prepost_1...
    (time_window, pre_visit1_time, iti_time, post_visit4_time);
%}

%[rm, rm_z, ctxid, visit_idx, half_idx, all_cell_idx_timefcn] = ALL_countwindow_prepost_1...
 %    (time_window, pre_visit1_time, iti_time, post_visit4_time, all_cell_idx);


%calculate population distance to same context type as a function of visit
%number
dists_out = dist_to_same_bw(time_window);
%dists_out = dist_to_same_bw_reduced(time_window, all_cell_idx_bw);
bins = 30;corr_to_same_bw(bins)

%anova text with dists out (above) calculated seperately for within and between
%
anovan_y = [dists_out_within{1}; dists_out_within{2}; dists_out_within{3}; dists_out_between{1}; dists_out_between{2}; dists_out_between{3}];
g1 = [ones(size([dists_out_within{1}; dists_out_within{2}; dists_out_within{3}])); repmat(2, size([dists_out_between{1}; dists_out_between{2}; dists_out_between{3}]))];
g2 = [ones(size(dists_out_within{1})); repmat(2, size(dists_out_within{2})); repmat(3, size(dists_out_within{3})); ones(size(dists_out_between{1})); repmat(2, size(dists_out_between{2})); repmat(3, size(dists_out_between{3}))];
[p,tbl,stats] = anovan(anovan_y, {g1, g2}, 'model', 'interaction');
%}



%local time classification
%
%first two trials, second two trials
clearvars data_hold all_visidx
%data
%all_data_hold = rm_z(ismember(visit_idx(:,4), [1 2 3 4]),:);
data_hold{1} = rm_z(visit_idx(:,4)==1, :); 
data_hold{2} = rm_z(visit_idx(:,4)==2, :);
data_hold{3} = rm_z(visit_idx(:,4)==3, :); 
data_hold{4} = rm_z(visit_idx(:,4)==4, :);

%visits
all_visidx_1 = visit_idx(visit_idx(:,4)==1, 5); 
    all_visidx{1} = all_visidx_1 - repmat(min(all_visidx_1), size(all_visidx_1));
all_visidx_2 = visit_idx(visit_idx(:,4)==2, 5);
    all_visidx{2} = all_visidx_2 - repmat(min(all_visidx_2), size(all_visidx_2));
all_visidx_3 = visit_idx(visit_idx(:,4)==3, 5);
    all_visidx{3} = all_visidx_3 - repmat(min(all_visidx_3), size(all_visidx_3));
all_visidx_4 = visit_idx(visit_idx(:,4)==4, 5);
    all_visidx{4} = all_visidx_4 - repmat(min(all_visidx_4), size(all_visidx_4));

train_1 = [1 3];
test_1 = [2 4];
data_train_ia_1 = [data_hold{train_1(1)}; data_hold{train_1(2)}];
data_test_ia_1 = [data_hold{test_1(1)}; data_hold{test_1(2)}];
class_train_ia_1 = [all_visidx{train_1(1)}; all_visidx{train_1(2)}];
class_test_ia_1 = [all_visidx{test_1(1)}; all_visidx{test_1(2)}];

train_2 = [2 4];
test_2 = [1 3];
data_train_ia_2 = [data_hold{train_2(1)}; data_hold{train_2(2)}];
data_test_ia_2 = [data_hold{test_2(1)}; data_hold{test_2(2)}];
class_train_ia_2 = [all_visidx{train_2(1)}; all_visidx{train_2(2)}];
class_test_ia_2 = [all_visidx{test_2(1)}; all_visidx{test_2(2)}];

[p_correct_localtime_1, assignments_localtime_1, distances_localtime_1, class_test_localtime_1] = md_classify_traintest(data_train_ia_1, class_train_ia_1, data_test_ia_1, class_test_ia_1);
[p_correct_localtime_2, assignments_localtime_2, distances_localtime_2, class_test_localtime_2] = md_classify_traintest(data_train_ia_2, class_train_ia_2, data_test_ia_2, class_test_ia_2);
