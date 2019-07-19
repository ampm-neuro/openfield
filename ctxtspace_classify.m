function [p_correct, assignments, classinput_ts_idx, correct_spaceerror, incorrect_spaceerror, chance_spaceerror] = ctxtspace_classify(spk_cts, ctxt_idx, posx, posy, bins)
%takes spikecounts from timewindows, context ids, and positions, and outputs classifications
%and success rate from classifying tws into spatial bins

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
%try
            mtx_lin_idx = sub2ind([bins,bins], bin_idx_x_hold, bin_idx_y_hold);
            mtx_lin_idx_hold = nan(size(bin_idx_x)); mtx_lin_idx_hold(rng_idx) = mtx_lin_idx;
            mtx_lin_idx = mtx_lin_idx_hold;
%catch
%figure; plot(bin_idx_x_hold); title x
%figure; plot(bin_idx_y_hold); title y
%error('duh')
%end

            %context correction (number of bins times context visit minus 1)
            ctxt_correction = (find(rctxts == ictxt) -1) * bins^2;

            %load timespace idex by context and cluster
            timespace_idx(ctxt_idx==ictxt, iclust) = mtx_lin_idx + repmat(ctxt_correction, size(mtx_lin_idx));

        end
    end
    

%find minimum number of visits to each bin for each cell
for iclust = 1:size(spk_cts,2)
    bin_cts(:,iclust) = histcounts(timespace_idx(:,iclust), 1:(4*bins^2+1))';
end

%figure; imagesc(bin_cts)
%figure; histogram(bin_cts(:))

%remove cells with insufficient visits to at least one bin
%
min_visits = 40;
number_of_dropped_cells = sum(min(bin_cts)<min_visits)
spk_cts(:, min(bin_cts)<min_visits) = [];
timespace_idx(:, min(bin_cts)<min_visits) = [];
posx(:, min(bin_cts)<min_visits) = [];
posy(:, min(bin_cts)<min_visits) = [];
bin_cts(:, min(bin_cts)<min_visits) = [];
%}

classinput_ts_idx = [];
classinput_s_cts = [];

min_visits_hold = nan(((bins^2) * 4),1);

%DO THIS

%first min_visits to each bin
for ib = 1:((bins^2) * 4)
    
    %smallest number of visits to that bin over all cells (sessions)
    %min_visits = min(bin_cts(ib,:));
    min_visits_hold(ib) = min_visits;
    
    %add bin to new timespace index
    classinput_ts_idx = [classinput_ts_idx; repmat(ib, min_visits, 1)];
    
    %iterate though cells taking the first 'min_visits' number of visits
    current_bin_scs = nan(min_visits, size(spk_cts,2));
    for ic = 1:size(spk_cts,2)
        
        %current cell
        cbsc_hold = spk_cts(:,ic);
        
        %visits to that bin
        cbsc_hold = cbsc_hold(timespace_idx(:,ic)==ib);
        
        %load first 'min_visits' number of visits
        %current_bin_scs(:,ic) = cbsc_hold(1:min_visits);
        current_bin_scs(:,ic) = cbsc_hold(round(linspace(1,length(cbsc_hold),min_visits)));
    end
    
    classinput_s_cts = [classinput_s_cts; current_bin_scs];
    
end

%figure; histogram(min_visits_hold)


%classify
%set 4 contexts to 2
bin_ranges = [1 : bins^2; ...
    bins^2+1 : 2*bins^2; ...
    2*bins^2+1 : 3*bins^2; ...
    3*bins^2+1 : 4*bins^2];

classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:)))...
    - repmat(bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(2,:)))));%all black
classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:)))...
    - repmat(bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(3,:)))));%all white 1
classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:))) = classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:)))...
    - repmat(2*bins^2, size(classinput_ts_idx(ismember(classinput_ts_idx, bin_ranges(4,:)))));%all white 2
%mdclass
[p_correct, assignments] = md_classify(classinput_s_cts, classinput_ts_idx);



%positions in space ignoring time and context
%{
corrected_class = classinput_ts_idx;
corrected_assign = assignments;
for ctxt = 4:-1:2
    
    %subtract items in class and assign to reduce them all to
    %"first-context" linear indices
    context_correction = (ctxt -1) * bins^2;
    corrected_class(corrected_class > context_correction) = corrected_class(corrected_class > context_correction) - context_correction;
    corrected_assign(corrected_assign > context_correction) = corrected_assign(corrected_assign > context_correction) - context_correction;

end
%}


%distances during correct classifications
correct_spaceerror = [];
for ctxt = 1:2
    
    %all bins in this context
    class_context = [1:bins^2; bins^2+1:bins^2*2];
    
    %index bins class in context
    ccidx = ismember(classinput_ts_idx, class_context(ctxt,:));
    
    %index bins assignments in context
    caidx = ismember(assignments, class_context(ctxt,:));
    
    %cclass (context specific classifications)
    
    if length(classinput_ts_idx) > length(ccidx & caidx)
        size(classinput_ts_idx)
        size(ccidx)
        size(caidx)
        error('wtf')
    end
    
    cclass = classinput_ts_idx(ccidx & caidx);
    %ccassign
    ccassign = assignments(ccidx & caidx);
    
    
    %set all bins to first context for absolute position
    if any(cclass > bins^2)
        cclass = cclass - repmat(bins^2, size(cclass));
    end
    if any(ccassign > bins^2)
        ccassign = ccassign - repmat(bins^2, size(ccassign));
    end
    
    %change linear indices into row,column indices
    [I_actual, J_actual] = ind2sub([bins, bins], cclass);
    [I_decode, J_decode] = ind2sub([bins, bins], ccassign);
    
    
    

    %calculate distance between decoded and actual at each time_window
    correct_spaceerror = [correct_spaceerror; diag(dist([I_decode,J_decode],[I_actual,J_actual]'))];
    
end

%distances during incorrect classifications
incorrect_spaceerror = [];
for ctxt = 1:2
    
    %all bins in this context
    class_context = [1:bins^2; bins^2+1:bins^2*2];
    
    %index bins class in context
    ccidx = ismember(classinput_ts_idx, class_context(ctxt,:));
    
    %index bins assignments in context
    caidx = ismember(assignments, class_context(ctxt,:));
    
    %cclass (context specific classifications)
    cclass = classinput_ts_idx(ccidx & ~caidx);
    %ccassign
    ccassign = assignments(ccidx & ~caidx);
    
    %set all bins to first context for absolute position
    if any(cclass > bins^2)
        cclass = cclass - repmat(bins^2, size(cclass));
    end
    if any(ccassign > bins^2)
        ccassign = ccassign - repmat(bins^2, size(ccassign));
    end
    
    
    %change linear indices into row,column indices
    [I_actual, J_actual] = ind2sub([bins, bins], cclass);
    [I_decode, J_decode] = ind2sub([bins, bins], ccassign);

    %calculate distance between decoded and actual at each time_window
    incorrect_spaceerror = [incorrect_spaceerror; diag(dist([I_decode,J_decode],[I_actual,J_actual]'))];
    
end    

%
%distances during chance classifications

num_shufs = 10;
chance_spaceerror = nan(length(classinput_ts_idx), num_shufs);
chance_histcounts = nan(num_shufs, ceil(pdist([1 1; bins bins])));
for shuf = 1:num_shufs
    for ctxt = 1:2

        %all bins in this context
        class_context = [1:bins^2; bins^2+1:bins^2*2];

        %index bins class in context
        ccidx = ismember(classinput_ts_idx, class_context(ctxt,:));

        %index bins assignments in context
        %caidx = ismember(assignments, class_context(ctxt,:));

        %cclass (context specific classifications)
        cclass = classinput_ts_idx(ccidx);
        cclass_shuf = cclass(randperm(length(cclass)));
        %ccassign
        %ccassign = assignments(ccidx);
        %ccassign = ccassign(randperm(length(ccassign)));

        %set all bins to first context for absolute position
        if any(cclass > bins^2)
            cclass = cclass - repmat(bins^2, size(cclass));
        end
        %{
        if any(ccassign > bins^2)
            ccassign = ccassign - repmat(bins^2, size(ccassign));
        end
        %}


        %change linear indices into row,column indices
        [I_actual, J_actual] = ind2sub([bins, bins], cclass);
        [I_decode, J_decode] = ind2sub([bins, bins], cclass_shuf);

        %calculate distance between decoded and actual at each time_window
        if ctxt == 1
            idx_lo = 1; idx_hi = sum(ccidx);
        else
            idx_lo = idx_hi+1; idx_hi = length(classinput_ts_idx);
        end
        chance_spaceerror(idx_lo:idx_hi, shuf) = diag(dist([I_decode,J_decode],[I_actual,J_actual]'));
    end 
    %load histcounts for shuffle statistic and plots
    chance_histcounts(shuf, :) = histcounts(chance_spaceerror(:, shuf), 0:ceil(pdist([1 1; bins bins])), 'normalization', 'probability');
end

%sort
chance_histcounts = sort(chance_histcounts);
%calculate shuffle lines
shuf_mean = mean(chance_histcounts, 1);
shuf_hi = chance_histcounts(ceil(num_shufs*.975), :);
shuf_lo = chance_histcounts(ceil(num_shufs*.025), :);

%plot
figure; hold on; 
histogram(correct_spaceerror.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')
histogram(incorrect_spaceerror.*(1/bins), (0:ceil(pdist([1 1; bins bins]))).*(1/bins), 'normalization', 'probability')
%shuf range plot
xticks_shuf = (.5:ceil(pdist([1 1; bins bins]))-.5).*(1/bins);
plot(xticks_shuf, shuf_mean, 'k-', 'linewidth', 3)
plot(xticks_shuf, shuf_hi, 'k-', 'linewidth', 2)
plot(xticks_shuf, shuf_lo, 'k-', 'linewidth', 2)

%details
set(gca,'TickLength',[0, 0]); box off
xlabel('Distance (m)')







































end