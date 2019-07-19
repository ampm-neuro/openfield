function [all_within, all_between, all_shuffle] = ctxtspace_classify_crossbw(spk_cts, ctxt_idx, posx, posy, bins, min_visits)
%takes spikecounts from timewindows, context ids, and positions, and outputs classifications
%and success rate from classifying tws into spatial bins
%outputs are matrices of cells with 
%columns 
        %1 percent correct classifications;
        %2 spatial bin classification assignments;
        %3 distance between classifications and actual locations;
        %4 actual locations;
%and rows corresponding to the training and test groups (e.g., b1b2; see below)




%format ia to match bw
for iac = 5:8
    ctxt_idx(ctxt_idx==iac) = iac-4;
end

%relevant contexts
rctxts = unique(ctxt_idx(ctxt_idx<9))';

%cut timewindows from irrelevant contexts
rctxts_idx = ismember(ctxt_idx, rctxts);
spk_cts = spk_cts(rctxts_idx, :);
posx = posx(rctxts_idx, :);
posy = posy(rctxts_idx, :);
ctxt_idx = ctxt_idx(rctxts_idx, :);

%preallocate
timespace_idx = nan(size(spk_cts));
bin_cts = nan((bins^2) * 4, size(spk_cts,2));

%find spatial bin classifications
%0 and 1 are hard coded min and max for both x and y (see trials_pcia)
bin_edges = linspace(0, 1, bins+1);

%error prevention
posx(posx<0) = 0;
posx(posx>1) = 1;
posy(posy<0) = 0;
posy(posy>1) = 1;

%timespace bin index
%
    %for each context
    for ictxt = rctxts

        %for each cell (shorthand for session -> could be improved)
        for iclust = 1:size(spk_cts,2)

            %find bin indices for each x and... 
            posx_hold = posx(ctxt_idx==ictxt, iclust);
            [~,~, bin_idx_x] = histcounts(posx_hold, bin_edges);
            bin_idx_x(isnan(posx_hold)) = nan;
            %...y pos
            posy_hold = posy(ctxt_idx==ictxt, iclust);
            [~,~, bin_idx_y] = histcounts(posy_hold, bin_edges);
            bin_idx_y(isnan(posx_hold)) = nan;
            
            %exclude nans
            rng_idx = bin_idx_x>0 & bin_idx_x<=bins & bin_idx_y>0 & bin_idx_y<=bins;
            bin_idx_x_hold = bin_idx_x(rng_idx);
            bin_idx_y_hold = bin_idx_y(rng_idx);    

            %find matrix index (translate i,j to col vect number)
            mtx_lin_idx = sub2ind([bins,bins], bin_idx_x_hold, bin_idx_y_hold);
            mtx_lin_idx_hold = nan(size(bin_idx_x)); mtx_lin_idx_hold(rng_idx) = mtx_lin_idx;
            mtx_lin_idx = mtx_lin_idx_hold;

            %context correction (number of bins times context visit minus 1)
            ctxt_correction = (find(rctxts == ictxt) -1) * bins^2;

            %load timespace index by context and cluster
            timespace_idx(ctxt_idx==ictxt, iclust) = mtx_lin_idx + repmat(ctxt_correction, size(mtx_lin_idx));

        end
    end
    

%find minimum number of visits to each bin for each cell
for iclust = 1:size(spk_cts,2)
    bin_cts(:,iclust) = histcounts(timespace_idx(:,iclust), 1:(4*bins^2+1))';
end

%remove cells with insufficient visits to at least one bin
%{
%min_visits = 30;
number_of_dropped_cells = sum(min(bin_cts)<min_visits)
spk_cts(:, min(bin_cts)<min_visits) = [];
timespace_idx(:, min(bin_cts)<min_visits) = [];
posx(:, min(bin_cts)<min_visits) = [];
posy(:, min(bin_cts)<min_visits) = [];
%black_tms(:, min(bin_cts)<min_visits) = [];
%white_tms(:, min(bin_cts)<min_visits) = [];
bin_cts(:, min(bin_cts)<min_visits) = [];
%}

classinput_ts_idx = [];
classinput_s_cts = [];
min_visits_hold = nan(((bins^2) * 4),1);

%min_visits to each bin
ib_abs = repmat(1:bins^2, 1, 4);
for ib = 1:((bins^2) * 4)
    
    %smallest number of visits to that bin over all cells (sessions)
    %min_visits = min(bin_cts(ib,:));
    min_visits_hold(ib) = min_visits;
    
    %add bin to new timespace index
    classinput_ts_idx = [classinput_ts_idx; repmat(ib, min_visits, 1)]; %full timespace bins
    %classinput_ts_idx = [classinput_ts_idx; repmat(ib_abs(ib), min_visits, 1)]; %absolute position bins
    
    
    %iterate though cells taking the first 'min_visits' number of visits
    current_bin_scs = nan(min_visits, size(spk_cts,2)); %columns are rates during each visit for a particular cell, rows are vectors
    for ic = 1:size(spk_cts,2)
        
        %current cell
        cbsc_hold = spk_cts(:,ic);
        
        %visits to that bin
        cbsc_hold = cbsc_hold(timespace_idx(:,ic)==ib);
        
        %load 'min_visits' number of visits
        %current_bin_scs(:,ic) = cbsc_hold(1:min_visits);
        %NOT FIRST VISITS, BUT SPACED THROUGHOUT SESSION
        current_bin_scs(:,ic) = cbsc_hold(round(linspace(1,length(cbsc_hold),min_visits)));
        
    end
    
    classinput_s_cts = [classinput_s_cts; current_bin_scs]; %stacking for each spatial bin, rows now indicate visits with spatial bins chunks
end

%screen for nans at each time window
for ntw = 1:size(classinput_s_cts,1)
    if any(isnan(classinput_s_cts(ntw,:)))
        classinput_s_cts(:, isnan(classinput_s_cts(ntw,:))) = [];
        warning([num2str(sum(isnan(classinput_s_cts(ntw,:)))) ' cells dropped bc nan'])
    end
end



%classify
bin_ranges = [1 : bins^2; ...
    bins^2+1 : 2*bins^2; ...
    2*bins^2+1 : 3*bins^2; ...
    3*bins^2+1 : 4*bins^2];

%set 4 contexts to 2
%{
classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:)))...
    - repmat(bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:)))));%all black
classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:)))...
    - repmat(bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:)))));%all white 1
%classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:)))...
%    - repmat(2*bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:)))));%all white 2
%}
%mdclass

%bins ranges for each context
%{
ranges = [ 1 : 2*min_visits*bins^2;...
    2*min_visits*bins^2+1 : 4*min_visits*bins^2 ];
black_ranges = ranges(1,:);
white_ranges = ranges(2,:);

%training means
%
tms = nan(length(unique(classinput_ts_idx(1 : 2*min_visits*bins^2))), size(classinput_s_cts(1 : 2*min_visits*bins^2, :),2), 2);
for ctx = 1:2
    rng = ranges(ctx,:);
    
    class_vect = classinput_ts_idx(rng);
    spikes_mtx = classinput_s_cts(rng, :);
    classes = unique(class_vect)';
    
    for im = classes 
        tms(classes==im, :, ctx) = mean(spikes_mtx(class_vect==im, :));
    end
end
%}

%mdclass, correct training set

%black contexts
b1_sc = classinput_s_cts(ismember(classinput_ts_idx, bin_ranges(1,:)), :);
b1_class = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(1,:)));
b2_sc = classinput_s_cts(ismember(classinput_ts_idx, bin_ranges(2,:)), :);
b2_class = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:))); 
    b2_class = b2_class - repmat(min(bin_ranges(2,:))-1, size(b2_class));
w1_sc = classinput_s_cts(ismember(classinput_ts_idx, bin_ranges(3,:)), :);
w1_class = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:))); 
    w1_class = w1_class - repmat(min(bin_ranges(3,:))-1, size(w1_class));
w2_sc = classinput_s_cts(ismember(classinput_ts_idx, bin_ranges(4,:)), :);
w2_class = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:))); 
    w2_class = w2_class - repmat(min(bin_ranges(4,:))-1, size(w2_class));    
    
%WITHIN CONTEXT
%
    data_train = b1_sc;
    class_train = b1_class;
    data_test = b2_sc;
    class_test = b2_class;
    [p_correct_b1b2, assignments_b1b2, ~, class_test_b1b2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b1b2 = assign_spaceerror(class_test_b1b2, assignments_b1b2, bins); %internal function
        %load output    
        all_within{1,1} = p_correct_b1b2;
        all_within{1,2} = assignments_b1b2;
        all_within{1,3} = spaceerror_b1b2;
        all_within{1,4} = class_test_b1b2;

    data_train = b2_sc;
    class_train = b2_class;
    data_test = b1_sc;
    class_test = b1_class;
    [p_correct_b2b1, assignments_b2b1, ~, class_test_b2b1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b2b1 = assign_spaceerror(class_test_b2b1, assignments_b2b1, bins); %internal function
        %load output    
        all_within{2,1} = p_correct_b2b1;
        all_within{2,2} = assignments_b2b1;
        all_within{2,3} = spaceerror_b2b1;
        all_within{2,4} = class_test_b2b1;

    %white contexts
    data_train = w1_sc;
    class_train = w1_class;
    data_test = w2_sc;
    class_test = w2_class;
    [p_correct_w1w2, assignments_w1w2, ~, class_test_w1w2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w1w2 = assign_spaceerror(class_test_w1w2, assignments_w1w2, bins); %internal function
        %load output    
        all_within{3,1} = p_correct_w1w2;
        all_within{3,2} = assignments_w1w2;
        all_within{3,3} = spaceerror_w1w2;
        all_within{3,4} = class_test_w1w2;

    data_train = w2_sc;
    class_train = w2_class;
    data_test = w1_sc;
    class_test = w1_class;
    [p_correct_w2w1, assignments_w2w1, ~, class_test_w2w1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w2w1 = assign_spaceerror(class_test_w2w1, assignments_w2w1, bins); %internal function
        %load output    
        all_within{4,1} = p_correct_w2w1;
        all_within{4,2} = assignments_w2w1;
        all_within{4,3} = spaceerror_w2w1;
        all_within{4,4} = class_test_w2w1;
        
%BETWEEN CONTEXT
%
    data_train = b1_sc;
    class_train = b1_class;
    data_test = w2_sc;
    class_test = w2_class;
    [p_correct_b1w2, assignments_b1w2, ~, class_test_b1w2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b1w2 = assign_spaceerror(class_test_b1w2, assignments_b1w2, bins); %internal function
        %load output    
        all_between{1,1} = p_correct_b1w2;
        all_between{1,2} = assignments_b1w2;
        all_between{1,3} = spaceerror_b1w2;
        all_between{1,4} = class_test_b1w2;

    data_train = w2_sc;
    class_train = w2_class;
    data_test = b1_sc;
    class_test = b1_class;
    [p_correct_w2b1, assignments_w2b1, ~, class_test_w2b1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w2b1 = assign_spaceerror(class_test_w2b1, assignments_w2b1, bins); %internal function
        %load output    
        all_between{2,1} = p_correct_w2b1;
        all_between{2,2} = assignments_w2b1;
        all_between{2,3} = spaceerror_w2b1;
        all_between{2,4} = class_test_w2b1;

    %white contexts
    data_train = w1_sc;
    class_train = w1_class;
    data_test = b2_sc;
    class_test = b2_class;
    [p_correct_w1b2, assignments_w1b2, ~, class_test_w1b2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w1b2 = assign_spaceerror(class_test_w1b2, assignments_w1b2, bins); %internal function
        %load output    
        all_between{3,1} = p_correct_w1b2;
        all_between{3,2} = assignments_w1b2;
        all_between{3,3} = spaceerror_w1b2;
        all_between{3,4} = class_test_w1b2;

    data_train = w2_sc;
    class_train = w2_class;
    data_test = b1_sc;
    class_test = b1_class;
    [p_correct_b2w1, assignments_b2w1, ~, class_test_b2w1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b2w1 = assign_spaceerror(class_test_b2w1, assignments_b2w1, bins); %internal function
        %load output    
        all_between{4,1} = p_correct_b2w1;
        all_between{4,2} = assignments_b2w1;
        all_between{4,3} = spaceerror_b2w1;
        all_between{4,4} = class_test_b2w1;
        
        
%SHUFFLE
%
    data_train = b1_sc;
    class_train = b1_class(randperm(length(b1_class)));
    data_test = b2_sc;
    class_test = b2_class;
    [p_correct_b1b2, assignments_b1b2, ~, class_test_b1b2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b1b2 = assign_spaceerror(class_test_b1b2, assignments_b1b2, bins); %internal function
        %load output    
        all_shuffle{1,1} = p_correct_b1b2;
        all_shuffle{1,2} = assignments_b1b2;
        all_shuffle{1,3} = spaceerror_b1b2;
        all_shuffle{1,4} = class_test_b1b2;

    data_train = b2_sc;
    class_train = b2_class(randperm(length(b2_class)));
    data_test = b1_sc;
    class_test = b1_class;
    [p_correct_b2b1, assignments_b2b1, ~, class_test_b2b1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_b2b1 = assign_spaceerror(class_test_b2b1, assignments_b2b1, bins); %internal function
        %load output    
        all_shuffle{2,1} = p_correct_b2b1;
        all_shuffle{2,2} = assignments_b2b1;
        all_shuffle{2,3} = spaceerror_b2b1;
        all_shuffle{2,4} = class_test_b2b1;

    %white contexts
    data_train = w1_sc;
    class_train = w1_class(randperm(length(w1_class)));
    data_test = w2_sc;
    class_test = w2_class;
    [p_correct_w1w2, assignments_w1w2, ~, class_test_w1w2] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w1w2 = assign_spaceerror(class_test_w1w2, assignments_w1w2, bins); %internal function
        %load output    
        all_shuffle{3,1} = p_correct_w1w2;
        all_shuffle{3,2} = assignments_w1w2;
        all_shuffle{3,3} = spaceerror_w1w2;
        all_shuffle{3,4} = class_test_w1w2;

    data_train = w2_sc;
    class_train = w2_class(randperm(length(w2_class)));
    data_test = w1_sc;
    class_test = w1_class;
    [p_correct_w2w1, assignments_w2w1, ~, class_test_w2w1] = md_classify_traintest(data_train, class_train, data_test, class_test);
    spaceerror_w2w1 = assign_spaceerror(class_test_w2w1, assignments_w2w1, bins); %internal function
        %load output    
        all_shuffle{4,1} = p_correct_w2w1;
        all_shuffle{4,2} = assignments_w2w1;
        all_shuffle{4,3} = spaceerror_w2w1;
        all_shuffle{4,4} = class_test_w2w1;
        
        
   

%{
[p_correct_bb, assignments_bb] = md_classify(classinput_s_cts(black_ranges, :), classinput_ts_idx(black_ranges));
[p_correct_ww, assignments_ww] = md_classify(classinput_s_cts(white_ranges, :), classinput_ts_idx(white_ranges));
%[p_correct_bb, assignments_bb] = md_classify_inputmeans(classinput_s_cts(black_ranges, :), classinput_ts_idx(black_ranges), tms(:, :, 1));
%[p_correct_ww, assignments_ww] = md_classify_inputmeans(classinput_s_cts(white_ranges, :), classinput_ts_idx(white_ranges), tms(:, :, 2));
p_correct_within = mean([p_correct_bb p_correct_ww]);
assignments_within = [assignments_bb; assignments_ww];
spaceerror_within = assign_spaceerror(classinput_ts_idx, assignments_within, bins); %internal function

%mdclass, incorrect training set 
[p_correct_bw, assignment_bw] = md_classify_inputmeans(classinput_s_cts(black_ranges, :), classinput_ts_idx(black_ranges), tms(:, :, 2));
[p_correct_wb, assignments_wb] = md_classify_inputmeans(classinput_s_cts(white_ranges, :), classinput_ts_idx(white_ranges), tms(:, :, 1));
p_correct_between = mean([p_correct_bw p_correct_wb]);
assignments_between = [assignment_bw; assignments_wb];
spaceerror_between = assign_spaceerror(classinput_ts_idx, assignments_between, bins); %internal function

%mdclass, chance training set 
classinput_ts_idx_r = classinput_ts_idx(randperm(length(classinput_ts_idx))); %randomize spatial bin labels for TWs
[p_correct_bw_r, assignment_bw_r] = md_classify_inputmeans(classinput_s_cts(black_ranges, :), classinput_ts_idx_r(black_ranges), tms(:, :, 2));
[p_correct_wb_r, assignments_wb_r] = md_classify_inputmeans(classinput_s_cts(white_ranges, :), classinput_ts_idx_r(white_ranges), tms(:, :, 1));
p_correct_between_r = mean([p_correct_bw_r p_correct_wb_r]);
assignments_between_r = [assignment_bw_r; assignments_wb_r];
spaceerror_between_r = assign_spaceerror(classinput_ts_idx_r, assignments_between_r, bins); %internal function

p_correct_withinbetween = [p_correct_within p_correct_between];
assignments_withinbetween = [assignments_within assignments_between];
spaceerror_withinbetween = [spaceerror_within spaceerror_between];
%}
    %distances during correct classifications
    function spaceerror = assign_spaceerror(cclass, ccassign, bins)
        %change linear indices into row,column indices
        [I_actual, J_actual] = ind2sub([bins, bins], cclass);
        [I_decode, J_decode] = ind2sub([bins, bins], ccassign);
        %calculate distance between decoded and actual at each time_window
        spaceerror = diag(dist([I_decode,J_decode],[I_actual,J_actual]'));
    end

%sort
%{
chance_histcounts = sort(chance_histcounts);
%calculate shuffle lines
shuf_mean = mean(chance_histcounts, 1);
shuf_hi = chance_histcounts(ceil(num_shufs*.975), :);
shuf_lo = chance_histcounts(ceil(num_shufs*.025), :);
%}

%plot
figure; hold on; 
spaceerror_within = [all_within{1,3}; all_within{2,3}; all_within{3,3}; all_within{4,3}];
spaceerror_between = [all_between{1,3}; all_between{2,3}; all_between{3,3}; all_between{4,3}];
spaceerror_shuffle = [all_shuffle{1,3}; all_shuffle{2,3}; all_shuffle{3,3}; all_shuffle{4,3}];

histogram(spaceerror_within.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')
histogram(spaceerror_between.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')
histogram(spaceerror_shuffle.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')
%histogram(spaceerror_between_r.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')


%shuf range plot
%xticks_shuf = (.5:ceil(pdist([1 1; bins bins]))-.5).*(1/bins);
%{
plot(xticks_shuf, shuf_mean, 'k-', 'linewidth', 3)
plot(xticks_shuf, shuf_hi, 'k-', 'linewidth', 2)
plot(xticks_shuf, shuf_lo, 'k-', 'linewidth', 2)
%}

%details
set(gca,'TickLength',[0, 0]); box off
xlabel('Distance (m)')
%}






































end