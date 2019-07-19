%pull out populations
load('C:\Users\ampm1\Desktop\oldmatlab\PCia\pop_ctx_9_27_simple.mat')
load('C:\Users\ampm1\Desktop\oldmatlab\PCia\250ms_with_hd.mat')
%fraction of interest
foi = 1/2;

%number of shuffles
shufs = 100;

all_distances = [];
% (remove columns)
posterior_for_vid = [];
bw_pop_cell_idx_hold = bw_pop_cell_idx;
spk_cts = all_counts_bw;
spk_cts = spk_cts(:, bw_pop_cell_idx_hold > 0);
posx = all_hd_bw(:,bw_pop_cell_idx_hold > 0,1);
bw_pop_cell_idx_hold = bw_pop_cell_idx_hold(bw_pop_cell_idx_hold > 0);

%cut timewindows from irrelevant contexts
% (remove rows)
rctxts = unique(context_bw(context_bw<9))'; %1:4
rctxts_idx = ismember(context_bw, rctxts);
spk_cts = spk_cts(rctxts_idx, :);
posx = posx(rctxts_idx, :);
ctxt_idx = context_bw(rctxts_idx, :);


%timespace bin index
%

%find spatial bin classifications
%0 and 360 are hard coded min and max for both x and y (see trials_pcia)
%bins = 6;
bins = 7; bins = bins^2;
bin_edges = linspace(0, 360, bins+1);

%error prevention
posx(posx<0) = 0;
posx(posx>360) = 360;


%context of interest
CoI = 3;

%for each context of interest
for ictxt = CoI

    %for each cell (shorthand for session -> could be improved)
    for iclust = 1:size(spk_cts,2)

        %find bin indices for each x and... 
        posx_hold = posx(ctxt_idx==ictxt, iclust);
        [~, ~, bin_idx_x] = histcounts(posx_hold, bin_edges);
        bin_idx_x(isnan(posx_hold)) = nan;

        %exclude nans
        rng_idx = bin_idx_x>0 & bin_idx_x<=bins;

        %load space index by context and cluster
        timespace_idx(ctxt_idx==ictxt, iclust) = bin_idx_x(rng_idx); 

    end
end

%preallocate context_matchup X population  
out_class_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
out_true_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
accuracy_wc = nan(1/foi, length(unique(bw_pop_cell_idx_hold)));


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
            included_pixels = 1:bins;

            %remove training pixles that weren't visited for at least 5sec.
            %
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
            %
            deletion_idx_test = ...
                ~ismember(group_test_wc, group_train_wc);
            sample_wc(deletion_idx_test, :) = [];
            group_test_wc(deletion_idx_test) = [];
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
            %figure; plot(class_wc)
            %hold on; plot(group_test_wc)

            %posterior_hold(:, included_pixels) = posterior_wc;
            posterior_wc = posterior_wc./sum(posterior_wc(:));

            %load
            out_class_wc{iwc, ipop} =  class_wc;
            out_true_wc{iwc, ipop} =  group_test_wc;
            accuracy_wc(iwc, ipop) = sum(class_wc==group_test_wc)/length(group_test_wc);

            %output posterior of largest population
            %{
            if size(spike_counts,2)==11
                posterior_for_vid = [posterior_for_vid; posterior_wc];
            end
            %}

        end

        %distances (in bins) between decoded positions and true positions
        all_class = [cell2mat(out_class_wc(1,ipop)); cell2mat(out_class_wc(2,ipop))];
        all_true = [cell2mat(out_true_wc(1,ipop)); cell2mat(out_true_wc(2,ipop))];

            for i = 1:length(all_true)
                all_distances = [all_distances; circ_distance(all_class(i), all_true(i), [1 bins])];
            end

        accuracy_distances(ipop) = sum(all_distances <= (100/bins)*sqrt(2) )/length(all_distances);

        %load shuffle catches
        abs_acc(ishuf, ipop) = mean(accuracy_wc(iwc, ipop));
        dist_acc(ishuf, ipop) = accuracy_distances(ipop);

    end
   
end


%average across contexts
shuffle_abs_acc = mean(abs_acc(2:end, :),2);
shuffle_abs_acc = sort(shuffle_abs_acc);
shuffle_dist_acc = mean(dist_acc(2:end, :),2);
shuffle_dist_acc = sort(shuffle_dist_acc);

%shuffle_range = [shuffle_dist_acc(ceil(.025*shufs)) mean(shuffle_dist_acc) shuffle_dist_acc(ceil(.975*shufs))];

shuffle_range = [shuffle_abs_acc(ceil(.025*shufs)) mean(shuffle_abs_acc) shuffle_abs_acc(ceil(.975*shufs))];

hold on;
for i = 1:3; hold on; plot(xlim, [shuffle_range(i) shuffle_range(i)], 'r-'); end







