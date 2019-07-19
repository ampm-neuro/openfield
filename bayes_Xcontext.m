%pull out populations
%bw_pop_cell_idx from ALL_cellcounts
% (remove columns)
bw_pop_cell_idx_hold = bw_pop_cell_idx;
spk_cts = all_counts_bw;
spk_cts = spk_cts(:, bw_pop_cell_idx_hold > 0);
posx = all_pos_bw(:,bw_pop_cell_idx_hold > 0,1);
posy = all_pos_bw(:,bw_pop_cell_idx_hold > 0,2);
bw_pop_cell_idx_hold = bw_pop_cell_idx_hold(bw_pop_cell_idx_hold > 0);

%cut timewindows from irrelevant contexts
% (remove rows)
rctxts = unique(context_bw(context_bw<9))';
rctxts_idx = ismember(context_bw, rctxts);
spk_cts = spk_cts(rctxts_idx, :);
posx = posx(rctxts_idx, :);
posy = posy(rctxts_idx, :);
ctxt_idx = context_bw(rctxts_idx, :);


%timespace bin index
%

%find spatial bin classifications
%0 and 1 are hard coded min and max for both x and y (see trials_pcia)
bins = 6;
bin_edges = linspace(0, 1, bins+1);

%error prevention
posx(posx<0) = 0;
posx(posx>1) = 1;
posy(posy<0) = 0;
posy(posy>1) = 1;

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
        %ctxt_correction = (find(rctxts == ictxt) -1) * bins^2;

        %load timespace index by context and cluster
        timespace_idx(ctxt_idx==ictxt, iclust) = mtx_lin_idx; %+ repmat(ctxt_correction, size(mtx_lin_idx));

    end
end


%within context matchups
wc = [1 2; 2 1; 3 4; 4 3];
%between context matchups
bc = [1 3; 3 1; 1 4; 4 1; 2 3; 3 2; 2 4; 4 2];


%preallocate context_matchup X population  
out_class_wc = cell(size(wc,1), length(unique(bw_pop_cell_idx_hold)));
out_true_wc = cell(size(wc,1), length(unique(bw_pop_cell_idx_hold)));
accuracy_wc = nan(size(wc,1), length(unique(bw_pop_cell_idx_hold)));

out_class_bc = cell(size(bc,1), length(unique(bw_pop_cell_idx_hold)));
out_true_bc = cell(size(bc,1), length(unique(bw_pop_cell_idx_hold)));
accuracy_bc = nan(size(wc,1), length(unique(bw_pop_cell_idx_hold)));

%bayesian decode each pop session
for ipop = unique(bw_pop_cell_idx_hold)'
    
    %for each within context matchup
    for iwc = 1:size(wc,1)
        
        %current matchup
        wc_c = wc(iwc,:);
    
        %training sample (context)
        training_wc = spk_cts(ctxt_idx == wc_c(1), bw_pop_cell_idx_hold == ipop);
        
        %true training class
        group_train_wc = timespace_idx(ctxt_idx == wc_c(1), bw_pop_cell_idx_hold == ipop);
        group_train_wc = group_train_wc(:,1); %columns are redundant
        
        %testing sample (context)
        sample_wc = spk_cts(ctxt_idx == wc_c(2), bw_pop_cell_idx_hold == ipop);
        
        %true test class
        group_test_wc = timespace_idx(ctxt_idx == wc_c(2), bw_pop_cell_idx_hold == ipop);
        group_test_wc = group_test_wc(:,1); %columns are redundant
        
        %decode
        [class_wc, posterior_wc] = bayesian_decode(sample_wc, training_wc, group_train_wc, .25);
    
        %load
        out_class_wc{iwc, ipop} =  class_wc;
        out_true_wc{iwc, ipop} =  group_test_wc;
        accuracy_wc(iwc, ipop) = sum(class_wc==group_test_wc);
        
    end
    
    
    
    %for each between context matchup
    for ibc = 1:size(bc,1)
        
        %current matchup
        bc_c = bc(ibc,:);
    
        %training sample (context)
        training_bc = spk_cts(ctxt_idx == bc_c(1), bw_pop_cell_idx_hold == ipop);
        
        %true training class
        group_train_bc = timespace_idx(ctxt_idx == bc_c(1), bw_pop_cell_idx_hold == ipop);
        group_train_bc = group_train_bc(:,1);
        
        %testing sample (context)
        sample_bc = spk_cts(ctxt_idx == bc_c(2), bw_pop_cell_idx_hold == ipop);
        
        %true test class
        group_test_bc = timespace_idx(ctxt_idx == bc_c(2), bw_pop_cell_idx_hold == ipop);
        group_test_bc = group_test_bc(:,1); %columns are redundant
        
        %decode
        [class_bc, posterior_bc] = bayesian_decode(sample_bc, training_bc, group_train_bc, .25);
        
        %load
        out_class_bc{ibc, ipop} =  class_bc;
        out_true_bc{ibc, ipop} =  group_test_bc;
        accuracy_bc(ibc, ipop) = sum(class_bc==group_test_bc);
    
    end
end

%plot
figure; hold on
cell_e{1} = mean(accuracy_wc)';
cell_e{2} =  mean(accuracy_bc)';
plot([ones(size(mean(accuracy_wc)))' (ones(size(mean(accuracy_bc))).*2)']',...
    [mean(accuracy_wc)' mean(accuracy_bc)']', '-', 'color', [.8 .8 .8])
errorbar_plot(cell_e)
ylim([0 inf])

%chance line
plot([0.5 2.5], [1/bins^2 1/bins^2].*length(out_class_wc{1,1}), 'k--')



