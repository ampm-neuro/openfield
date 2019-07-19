%pull out populations
%bw_pop_cell_idx from ALL_cellcounts
%or
load('C:\Users\ampm1\Desktop\oldmatlab\PCia\pop_ctx_9_27_simple.mat')
%fraction of interest
foi = 1/2;

iscores = [];

% (remove columns)
posterior_for_vid = [];
bw_pop_cell_idx_hold = bw_pop_cell_idx;
spk_cts = all_counts_bw;
spk_cts = spk_cts(:, bw_pop_cell_idx_hold > 0);
posx = all_pos_bw(:,bw_pop_cell_idx_hold > 0,1);
posy = all_pos_bw(:,bw_pop_cell_idx_hold > 0,2);
bw_pop_cell_idx_hold = bw_pop_cell_idx_hold(bw_pop_cell_idx_hold > 0);

%cut timewindows from irrelevant contexts
% (remove rows)
rctxts = unique(context_bw(context_bw<9))'; %1:4
rctxts_idx = ismember(context_bw, rctxts);
spk_cts = spk_cts(rctxts_idx, :);
posx = posx(rctxts_idx, :);
posy = posy(rctxts_idx, :);
ctxt_idx = context_bw(rctxts_idx, :);


%timespace bin index
%

%find spatial bin classifications
%0 and 1 are hard coded min and max for both x and y (see trials_pcia)
%bins = 6;
bins = 40;
bin_edges = linspace(0, 1, bins+1);

%error prevention
posx(posx<0) = 0;
posx(posx>1) = 1;
posy(posy<0) = 0;
posy(posy>1) = 1;


%context of interest
CoI = 3;

%for each context of interest
for ictxt = CoI

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

        %load space index by context and cluster
        timespace_idx(ctxt_idx==ictxt, iclust) = mtx_lin_idx; 

    end
end




%preallocate context_matchup X population  
out_class_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
out_true_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
accuracy_wc = nan(1/foi, length(unique(bw_pop_cell_idx_hold)));

%bayesian decode each pop session
for ipop = unique(bw_pop_cell_idx_hold)'
    
    %context spike counts
    spike_counts = spk_cts(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
    %context classifications
    classifications = timespace_idx(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
    classifications = classifications(:,1); %columns are redundant
    
    %fill cell for iteratations below
    %fraction of interest
    train_test = cell(1/foi, 2);
    for i = 1:1/foi
        foi_indx = floor(length(classifications)/(1/foi)*(i-1)+1):floor(length(classifications)/(1/foi)*i);
        train_test{i,1} = spike_counts(foi_indx, :);
        train_test{i,2} = classifications(foi_indx);
    end

    
    %for each half
    for iwc = 1:1/foi

        %training sample (context)
        training_wc = cat(1,train_test{setdiff(1:1/foi, iwc), 1});
        
        %true training class
        group_train_wc = cat(1,train_test{setdiff(1:1/foi, iwc), 2});
        %num_train_pix = length(unique(group_train_wc))

        %track excluded pixles
        included_pixels = unique(group_train_wc);
        
        %remove training pixles that weren't visited for at least 5sec.
        %{
        [a, b] = histc(group_train_wc, unique(group_train_wc));
        deletion_idx = ismember(b, find(a<(5*4)));
        group_train_wc(deletion_idx) = [];
        training_wc(deletion_idx, :) = [];
        included_pixels = intersect(included_pixels, unique(group_train_wc));
        %}
        
        %testing sample (context)
        sample_wc = train_test{iwc, 1};
        
        %true test class
        group_test_wc = train_test{iwc, 2}; 
        %num_test_pix = length(unique(group_test_wc))

        %remove test pixles that aren't in train set
        %{
        deletion_idx_test = ...
            ~ismember(group_test_wc, group_train_wc);
        sample_wc(deletion_idx_test, :) = [];
        group_test_wc(deletion_idx_test) = [];
        %}
        
        
        %reliability plots
        %
        if ipop==5
           
           f_x = nan(bins^2, size(training_wc,2));
           for pixel = unique(group_train_wc)'
                f_x(pixel, :) = mean(training_wc(group_train_wc==pixel, :))/0.25;
           end

           for icell = 1:size(training_wc,2)
            figure; 
            imagesc(inpaint_nans(smooth2a(reshape(f_x(:,icell), bins, bins), 3))); 
            title(['cell ' num2str(icell)]); 
            set(gca,'TickLength',[0, 0]);
            axis square
            axis off
            colorbar
            colormap jet
           end
           
        end
        %}
        
        
        
        
        %{
        testclass = unique(group_test_wc);
        trainclass = unique(group_train_wc);
        all_pxl = zeros(7,7);
        all_pxl(trainclass) = 1;
        all_pxl(testclass) = 2;
        figure; imagesc(all_pxl); colorbar; caxis([0 2]);
        %}
        
        %decode
        [class_wc, posterior_wc] = bayesian_decode(sample_wc, training_wc, group_train_wc, .25);
        
        posterior_hold = nan(size(posterior_wc, 1), bins^2);
        posterior_hold(:, included_pixels) = posterior_wc;

        
        %guess missing values
        %{
        for iframe = 1:size(posterior_hold,1)
           frame_hold = reshape(posterior_hold(iframe,:), bins, bins);
           frame_hold = smooth2a(frame_hold, 2);
           frame_hold = inpaint_nans(frame_hold);
           frame_hold = frame_hold(:);
           posterior_hold(iframe, :) = frame_hold;
           class_wc(iframe) = find(frame_hold == max(frame_hold), 1, 'first');
        end
        %}
        
        posterior_wc = posterior_hold./sum(posterior_hold(:));
        
        %load
        out_class_wc{iwc, ipop} =  class_wc;
        out_true_wc{iwc, ipop} =  group_test_wc;
        accuracy_wc(iwc, ipop) = sum(class_wc==group_test_wc)/length(group_test_wc);
        
        %output posterior of largest population
        
        if size(spike_counts,2)==11
            posterior_for_vid = [posterior_for_vid; posterior_wc];
        end
        
    end
    
    %distances between decoded positions and true positions
    all_class = [cell2mat(out_class_wc(1,ipop)); cell2mat(out_class_wc(2,ipop))];
        [all_class_i, all_class_j] = ind2sub([bins bins], all_class);
        all_class = [all_class_i all_class_j];
    all_true = [cell2mat(out_true_wc(1,ipop)); cell2mat(out_true_wc(2,ipop))];
        [all_true_i, all_true_j] = ind2sub([bins bins], all_true);
        all_true = [all_true_i all_true_j];
    all_distances = dist(all_class, all_true');
    all_distances = all_distances(diag_mask(size(all_class,1)));
    
    %set distances to cm
    all_distances = all_distances.*(100/bins);
    accuracy_distances(ipop) = sum(all_distances <= (100/bins)*sqrt(2) )/length(all_distances);
   
end


%plot
%
figure; hold on
cell_e{1} = mean(accuracy_wc)';
errorbar_plot(cell_e)
ylim([0 inf])
plot([0.5 2.5], [1/bins^2 1/bins^2], 'k--')
title absolute
%}
%{
figure; hold on
cell_e{1} = accuracy_distances';
errorbar_plot(cell_e)
ylim([0 inf])
title distances
%}





