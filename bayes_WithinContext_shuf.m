%pull out populations
%bw_pop_cell_idx from ALL_cellcounts
%or
%load('C:\Users\ampm1\Desktop\oldmatlab\PCia\pop_ctx_9_27_simple.mat')

shufs = 1000;

% (remove columns)
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
bins = 7;
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
out_class_wc = cell(2, length(unique(bw_pop_cell_idx_hold)));
out_true_wc = cell(2, length(unique(bw_pop_cell_idx_hold)));
accuracy_wc = nan(2, length(unique(bw_pop_cell_idx_hold)));


%preallocate shuffle catch
abs_acc = nan(shufs, length(unique(bw_pop_cell_idx_hold)));
dist_acc = nan(shufs, length(unique(bw_pop_cell_idx_hold)));

%bayesian decode each pop session
for ipop = unique(bw_pop_cell_idx_hold)'
    
    ipop
    
    %context spike counts
    spike_counts = spk_cts(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
    num_spikes = size(spike_counts,2);
    num_time_bins = size(spike_counts,1);
    
    %shuffles
    for ishuf = 1:shufs

        %shuffle cell's spike-counts between time bins
        for ispk = 1:num_spikes
            spike_counts(:,ispk) = spike_counts(randperm(num_time_bins),ispk);
        end
    
        %context classifications
        classifications = timespace_idx(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
        classifications = classifications(:,1); %columns are redundant

        %first half
        sc_first = spike_counts(1:floor(size(spike_counts,1)/2), :);
        class_first = classifications(1:floor(size(classifications,1)/2));

        %second half
        sc_second = spike_counts(ceil(size(spike_counts,1)/2):end, :);
        class_second = classifications(ceil(size(classifications,1)/2):end);

        %fill cell for iteratations below
        train_test{1,1} = sc_first; %first half spikes
        train_test{1,2} = class_first; %first half classifications
        train_test{2,1} = sc_second; %second half spikes
        train_test{2,2} = class_second; %second half classifications


        %for each half
        for iwc = 1:2

            %training sample (context)
            training_wc = train_test{iwc, 1};

            %true training class
            group_train_wc = train_test{iwc, 2};

            %testing sample (context)
            sample_wc = train_test{setdiff([1 2], iwc), 1};

            %true test class
            group_test_wc = train_test{setdiff([1 2], iwc), 2};

            %decode
            [class_wc, posterior_wc] = bayesian_decode(sample_wc, training_wc, group_train_wc, .25);

            %load
            out_class_wc{iwc, ipop} =  class_wc;
            out_true_wc{iwc, ipop} =  group_test_wc;
            accuracy_wc(iwc, ipop) = sum(class_wc==group_test_wc)/length(group_test_wc);

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
        
        %load shuffle catches
        abs_acc(ishuf, ipop) = mean(accuracy_wc(iwc, ipop));
        dist_acc(ishuf, ipop) = accuracy_distances(ipop);
        
   
    end
    
end

%average across contexts
shuffle_abs_acc = mean(abs_acc,2);
shuffle_abs_acc = sort(shuffle_abs_acc);
shuffle_dist_acc = mean(dist_acc,2);
shuffle_dist_acc = sort(shuffle_dist_acc);

shuffle_range = [shuffle_dist_acc(floor(.025*shufs)) mean(shuffle_dist_acc) shuffle_dist_acc(ceil(.975*shufs))];
hold on;
for i = 1:3; hold on; plot(xlim, [shuffle_range(i) shuffle_range(i)], 'k-'); end






