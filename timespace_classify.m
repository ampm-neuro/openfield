function [p_correct, assignments, classinput_ts_idx, spatial_error] = timespace_classify(spk_cts, ctxt_idx, posx, posy, bins, min_visits)
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

figure; histogram(bin_cts(:))

%remove cells with insufficient visits to at least one bin
%
%min_visits = 40;
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
[p_correct, assignments] = md_classify(classinput_s_cts, classinput_ts_idx);



%positions in space ignoring time and context
absolute_positions = [];
corrected_class = classinput_ts_idx;
corrected_assign = assignments;
for ctxt = 4:-1:2
    
    %subtract items in class and assign to reduce them all to
    %"first-context" linear indices
    context_correction = (ctxt -1) * bins^2;
    corrected_class(corrected_class > context_correction) = corrected_class(corrected_class > context_correction) - context_correction;
    corrected_assign(corrected_assign > context_correction) = corrected_assign(corrected_assign > context_correction) - context_correction;

end

%change linear indices into row,column indices
[I_actual,J_actual] = ind2sub([bins, bins], corrected_class);
[I_decode,J_decode] = ind2sub([bins, bins], corrected_assign);

%calculate distance between decoded and actual at each time_window
spatial_error = diag(dist([I_decode,J_decode],[I_actual,J_actual]'));
    





































end