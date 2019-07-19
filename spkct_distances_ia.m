function [within_comps, between_comps, a_bar, o_bar, ao_bar] = spkct_distances_ia(spike_counts, contexts, within_comparisons, between_comparisons, method1, method2, time_window)
% distances(spike_counts, contexts, comparisons, method1, method2)
%
% calculates distances between HD clouds of spikecounts computed with 
% ALL_countwindow
%
% within_comparisons should be a 1 column matrix
% between_comparisons should be a 2 column matrix
%
% method1 input should be 1 for means and 2 for clouds
% method2 input should be 1 for means and 2 for samples
%   thus:
%       1 and 1 calculate distance between means and means
%       1 and 2 calculate distance between means and samples
%       2 and 1 calculate distance between clouds and means
%       2 and 2 calculate distance between clouds and samples
%
% for ia run iteratively with just i i and then just a a comparisons in
% both within and between. Compare within output between iterations.
%
%
ao_bar = [];
a_bar = [];
o_bar = [];

%dimenionality correction
dc = sqrt(size(spike_counts,2));

%identify relevant contexts
unq_w = unique(within_comparisons); unq_w = unq_w(:);
unq_b = unique(between_comparisons); unq_b = unq_b(:);
rlvt_ctxts = unique([unq_w; unq_b]);

%constrain spike_counts to samples from relevant contexts
spike_counts = spike_counts(ismember(contexts, rlvt_ctxts), :);
contexts = contexts(ismember(contexts, rlvt_ctxts));

%distances to means
if method1 == 1

    %means to means (euclidean)
    if method2 == 1
        
        %calculate means
        ctxt_means = nan(length(rlvt_ctxts), size(spike_counts,2));
        for cm = 1:size(ctxt_means,1)
            ctxt_means(cm, :) = mean(spike_counts(contexts == rlvt_ctxts(cm),:));
        end
        
        %withins
        within_comps = nan(1, size(within_comparisons,1));
        for wc = 1: size(within_comps,1)
            current_contexts = within_comparisons(wc,:);
            within_comps(wc) = pdist([ctxt_means(rlvt_ctxts==current_contexts(1), :); ctxt_means(rlvt_ctxts==current_contexts(end), :)], 'euclidean')./dc;
        end

        %betweens
        between_comps = nan(size(between_comparisons,1));
        for bc = 1: size(between_comps,1)
            current_contexts = between_comparisons(bc,:);
            between_comps(bc) = pdist([ctxt_means(rlvt_ctxts==current_contexts(1), :); ctxt_means(rlvt_ctxts==current_contexts(end), :)], 'euclidean')./dc;
        end
        
        %figure within
        figure; hold on; 
        title withins
        xticks(1:size(within_comparisons,1))
        xlabel contexts
        ylabel distance
        for iw = 1:size(within_comparisons,1)
            bar(iw, within_comps(iw))
            xt=get(gca,'xticklabel');
            xt{iw}=num2str(within_comparisons(iw,:));
            set(gca,'xticklabel',xt)
        end
        %figure between
        figure; hold on; 
        title betweens
        xticks(1:size(between_comparisons,1))
        xlabel contexts
        ylabel distance
        for ib = 1:size(between_comparisons,1)
            bar(ib, between_comps(ib))            
            xt=get(gca,'xticklabel');
            xt{ib}=num2str(between_comparisons(ib,:));
            set(gca,'xticklabel',xt)
        end
        
        
        
        
        
        
    %SAMPLES TO MEANS (iterative euclidean) 
    %
    % Finds the euclid distance from every sample to each population mean
    % (recalculated sans that sample), and then averages the distances to the
    % "same" context (within) and the distances the the other contexts
    % (between).
    %
    % outputs one within and one between for each observation
    %
    elseif method2 == 2
        
        %preallocate
        within_comps = nan(length(rlvt_ctxts), size(spike_counts,1));
        between_comps = nan(length(rlvt_ctxts), size(spike_counts,1));
        
        %all samples
        allsamps = 1:size(spike_counts,1);
        
        %for each sample
        for samp = allsamps
            
            %spike_counts sans samp
            spike_counts_dropone = spike_counts(allsamps~=samp,:);
            context_ids_dropone = contexts(allsamps~=samp);
            
            %recalculate means for each samp (dropping test samp)
            ctxt_means_dropone = nan(length(rlvt_ctxts), size(spike_counts_dropone,2));
            for cm = 1:size(ctxt_means_dropone,1)
                ctxt_means_dropone(cm, :) = mean(spike_counts_dropone(context_ids_dropone == rlvt_ctxts(cm),:));
            end

            %calculate dist from samp to every context mean (also prod. extraneous)
            dist_to_means = pdist([spike_counts(samp,:); ctxt_means_dropone]); 
            dist_to_means = dist_to_means(1:size(ctxt_means_dropone,1))./dc;
            
            %within_comparisons for this sample's context
            cc_w = within_comparisons(logical(sum(within_comparisons == contexts(samp),2)),:);
            cc_w_comp_idx = rlvt_ctxts == cc_w;
            
            %between_comparisons for this sample's context
            cc_b = between_comparisons(logical(sum(between_comparisons == contexts(samp),2)),:);
            cc_b = cc_b(cc_b~=contexts(samp));
            cc_b_comp_idx = ismember(rlvt_ctxts, cc_b);
            
            %assign the within distance
            within_comps(cc_w_comp_idx, samp) = dist_to_means(cc_w_comp_idx);

            %assign the between distances
            between_comps(cc_b_comp_idx, samp) = dist_to_means(cc_b_comp_idx) -1;
     
        end
        
        %remove nans
        %within_comps = within_comps(~isnan(within_comps));
        %between_comps = between_comps(~isnan(between_comps));
        

        %figure
        figure; hold on; 
        title withins
        xticks(1:size(within_comparisons,1))
        xlabel contexts
        ylabel distance
        for iw = 1:size(within_comparisons,1)
            mean_hold = mean(within_comps(iw, contexts==rlvt_ctxts(iw)));
            bar(iw, mean_hold)
            errorbar(iw, mean_hold, std(within_comps(iw, contexts==rlvt_ctxts(iw)))./sqrt(length(within_comps(iw, contexts==rlvt_ctxts(iw)))), 'k.')
            xt=get(gca,'xticklabel');
            xt{iw}=num2str(within_comparisons(iw,:));
            set(gca,'xticklabel',xt)
        end

        figure; hold on; 
        title betweens
        xticks(1:size(between_comparisons,1))
        xlabel contexts
        ylabel distance
        for ib = 1:size(between_comparisons,1)
            
            bc = between_comparisons(ib,:); %compared contexts
            
            ctx_dists = [between_comps(rlvt_ctxts==bc(2), contexts==bc(1)) between_comps(rlvt_ctxts==bc(1), contexts==bc(2))]; %all relevant comparisons
            mean_hold = mean(ctx_dists); %mean
            
            if ismember(bc, [5 7; 6 8; 5 8; 6 7; 1 3; 2 4; 1 4; 2 3], 'rows')
                ao_bar = [ao_bar; ctx_dists'];
                continue
            elseif ismember(bc, [5 6])
                a_bar = [a_bar; ctx_dists'];
            elseif ismember(bc, [7 8])
                o_bar = [o_bar; ctx_dists'];
            end
            
            
            %plot
            bar(ib, mean_hold)
            errorbar(ib, mean_hold, std(ctx_dists)./sqrt(length(ctx_dists)), 'k.')
            
            xt=get(gca,'xticklabel');
            xt{ib}=num2str(bc);
            set(gca,'xticklabel',xt)
        end
        
        %custom bar
        %
        ib=2
        bar(ib+1, mean(ao_bar))
        errorbar(ib+1, mean(ao_bar), std(ao_bar)./sqrt(length(ao_bar)), 'k.')
        xt=get(gca,'xticklabel');
        xt{ib+1}='A&O';
        set(gca,'xticklabel',xt)
        %}
        
    end
            
        

    
%distances to clouds
elseif method1 == 2

    %MEANS TO CLOUDS (averaging forward and back mahalanobis)
    %
    % Finds the mahal distance from each mean to each other population cloud 
    % from the same (within) or different (between) context type. To find the
    % "within" distance, the mahal dist from each context mean to the
    % opposite cloud is averaged
    %
    % e.g. the mahal distance between c1 and c2 is the average of
    % c1mean to c2cloud   and    c2mean to c1cloud
    %
    % outputs one within_comp for each row in within_comparisons and one 
    % between_comp for each row in between_comparisons
    %
    if method2 == 1
        
        
        %calculate means and clouds
        ctxt_means = nan(length(rlvt_ctxts), size(spike_counts,2));
        ctxt_clouds = cell(length(rlvt_ctxts),1);
        for cm = 1:size(ctxt_means,1)
            ctxt_means(cm, :) = mean(spike_counts(contexts == rlvt_ctxts(cm),:));
            ctxt_clouds{cm} = spike_counts(contexts == rlvt_ctxts(cm),:);
        end
        
        %withins
        within_comps = nan(1, size(within_comparisons,1));
        for wc = 1: size(within_comps,1)
            
            current_contexts = within_comparisons(wc,:);
            
            mahal_forward_w = mahal(ctxt_means(rlvt_ctxts==current_contexts(1), :), ctxt_clouds{rlvt_ctxts==current_contexts(end)})./dc;
            mahal_backward_w = mahal(ctxt_means(rlvt_ctxts==current_contexts(end), :), ctxt_clouds{rlvt_ctxts==current_contexts(1)})./dc;
            
            %average forward and reverse 
            within_comps(wc) = mean([mahal_forward_w mahal_backward_w]);
        end

        %betweens
        between_comps = nan(size(between_comparisons,1));
        for bc = 1: size(between_comps,1)
            
            current_contexts = between_comparisons(bc,:);

            mahal_forward_b = mahal(ctxt_means(rlvt_ctxts==current_contexts(1), :), ctxt_clouds{rlvt_ctxts==current_contexts(end)})./dc;
            mahal_backward_b = mahal(ctxt_means(rlvt_ctxts==current_contexts(end), :), ctxt_clouds{rlvt_ctxts==current_contexts(1)})./dc;
            
            %average forward and reverse
            between_comps(bc) = mean([mahal_forward_b mahal_backward_b]);

        end
        
        
        %figure within
        figure; hold on; 
        title withins
        xticks(1:size(within_comparisons,1))
        xlabel contexts
        ylabel distance
        for iw = 1:size(within_comparisons,1)
            bar(iw, within_comps(iw))
            xt=get(gca,'xticklabel');
            xt{iw}=num2str(within_comparisons(iw,:));
            set(gca,'xticklabel',xt)
        end
        %figure between
        figure; hold on; 
        title betweens
        xticks(1:size(between_comparisons,1))
        xlabel contexts
        ylabel distance
        for ib = 1:size(between_comparisons,1)
            bar(ib, between_comps(ib))            
            xt=get(gca,'xticklabel');
            xt{ib}=num2str(between_comparisons(ib,:));
            set(gca,'xticklabel',xt)
        end
        
        
        
    %samples to clouds (iterative mahalanobis)
    %
    % Finds the mahal distance from every sample to each other population 
    % cloud.
    %
    % e.g., from observation 1 (member of context 1) to the clouds of
    % contexts 2 3 and 4.
    %
    % outputs one within and one between for each observation
    %
    elseif method2 == 2
        
        %preallocate
        within_comps = nan(length(rlvt_ctxts), size(spike_counts,1));
        between_comps = nan(length(rlvt_ctxts), size(spike_counts,1));
        
        %all samples
        allsamps = 1:size(spike_counts,1);
        
        %for each sample
        for samp = allsamps
            
            %spike_counts sans samp
            spike_counts_dropone = spike_counts(allsamps~=samp,:);
            contexts_dropone = contexts(allsamps~=samp);
            
            %preallocate
            cloud_dists = nan(length(rlvt_ctxts),1);
            
            %for each context
            for cds = 1:length(cloud_dists)
                
                %distance from samp to cloud
                cloud_dists(cds) = mahal(spike_counts(samp,:), spike_counts_dropone(contexts_dropone == rlvt_ctxts(cds),:))./dc;

            end
            
            %current withins and betweens
            %cc_w = within_comparisons(logical(sum(within_comparisons == contexts(samp),2)),:);
            %cc_b = within_comparisons(~logical(sum(within_comparisons == contexts(samp),2)),:);
            
            %within_comparisons for this sample's context
            cc_w = within_comparisons(logical(sum(within_comparisons == contexts(samp),2)),:);
            cc_w_comp_idx = rlvt_ctxts == cc_w;
            
            %between_comparisons for this sample's context
            cc_b = between_comparisons(logical(sum(between_comparisons == contexts(samp),2)),:);
            cc_b = cc_b(cc_b~=contexts(samp));
            cc_b_comp_idx = ismember(rlvt_ctxts, cc_b);

            %assign the within distance
            within_comps(cc_w_comp_idx, samp) = cloud_dists(cc_w_comp_idx);%mean(cloud_dists(cc_w));

            %assign the between distances
            between_comps(cc_b_comp_idx, samp) = cloud_dists(cc_b_comp_idx);%mean(cloud_dists(cc_b));
            
        end
        
        
        %figure
        figure; hold on; 
        title withins
        xticks(1:size(within_comparisons,1))
        xlabel contexts
        ylabel distance
        for iw = 1:size(within_comparisons,1)
            mean_hold = mean(within_comps(iw, contexts==rlvt_ctxts(iw)));
            bar(iw, mean_hold)
            errorbar(iw, mean_hold, std(within_comps(iw, contexts==rlvt_ctxts(iw)))./sqrt(length(within_comps(iw, contexts==rlvt_ctxts(iw)))), 'k.')
            xt=get(gca,'xticklabel');
            xt{iw}=num2str(within_comparisons(iw,:));
            set(gca,'xticklabel',xt)
        end

        figure; hold on; 
        title betweens
        xticks(1:size(between_comparisons,1))
        xlabel contexts
        ylabel distance
        for ib = 1:size(between_comparisons,1)
            
            bc = between_comparisons(ib,:); %compared contexts
            
            ctx_dists = [between_comps(rlvt_ctxts==bc(2), contexts==bc(1)) between_comps(rlvt_ctxts==bc(1), contexts==bc(2))]; %all relevant comparisons
            mean_hold = mean(ctx_dists); %mean
            
            if ismember(bc, [5 7; 6 8; 5 8; 6 7; 1 3; 2 4; 1 4; 2 3], 'rows')
                ao_bar = [ao_bar; ctx_dists'];
            end
            
            %plot
            bar(ib, mean_hold)
            errorbar(ib, mean_hold, std(ctx_dists)./sqrt(length(ctx_dists)), 'k.')
            xt=get(gca,'xticklabel');
            xt{ib}=num2str(bc);
            set(gca,'xticklabel',xt)
        end
        
        %custom bar
        %
        bar(ib+1, mean(ao_bar))
        errorbar(ib+1, mean(ao_bar), std(ao_bar)./sqrt(length(ao_bar)), 'k.')
        xt=get(gca,'xticklabel');
        xt{ib+1}='A&O';
        set(gca,'xticklabel',xt)
        %}
        
    end
end

%report
display(['within comparison is average of all distances between (1) context '...
    num2str(within_comparisons(1, 1)) ' and ' num2str(within_comparisons(1, end))...
    ' and (2) context ' num2str(within_comparisons(end, 1)) ' and ' num2str(within_comparisons(end, end))])

if size(between_comparisons,1) == 1
    display(['between comparison is average of all distances between context '...
        num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))])
elseif size(between_comparisons,1) == 2
    display(['between comparison is average of all distances between (1) context '...
        num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
        ' and (2) context ' num2str(between_comparisons(end, 1)) ' and ' num2str(between_comparisons(end, end))])
elseif size(between_comparisons,1) == 3
    display(['between comparison is average of all distances between (1) context '...
        num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
        ' and (2) context ' num2str(between_comparisons(2, 1)) ' and ' num2str(between_comparisons(2, end))...
        ' and (3) context ' num2str(between_comparisons(3, 1)) ' and ' num2str(between_comparisons(3, end))])
elseif size(between_comparisons,1) == 4
    display(['between comparison is average of all distances between (1) context '...
        num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
        ' and (2) context ' num2str(between_comparisons(2, 1)) ' and ' num2str(between_comparisons(2, end))...
        ' and (3) context ' num2str(between_comparisons(3, 1)) ' and ' num2str(between_comparisons(3, end))...
        ' and (4) context ' num2str(between_comparisons(4, 1)) ' and ' num2str(between_comparisons(4, end))])
elseif size(between_comparisons,1) == 5
display(['between comparison is average of all distances between (1) context '...
    num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
    ' and (2) context ' num2str(between_comparisons(2, 1)) ' and ' num2str(between_comparisons(2, end))...
    ' and (3) context ' num2str(between_comparisons(3, 1)) ' and ' num2str(between_comparisons(3, end))...
    ' and (4) context ' num2str(between_comparisons(4, 1)) ' and ' num2str(between_comparisons(4, end))...
    ' and (5) context ' num2str(between_comparisons(5, 1)) ' and ' num2str(between_comparisons(5, end))])
elseif size(between_comparisons,1) == 6
display(['between comparison is average of all distances between (1) context '...
    num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
    ' and (2) context ' num2str(between_comparisons(2, 1)) ' and ' num2str(between_comparisons(2, end))...
    ' and (3) context ' num2str(between_comparisons(3, 1)) ' and ' num2str(between_comparisons(3, end))...
    ' and (4) context ' num2str(between_comparisons(4, 1)) ' and ' num2str(between_comparisons(4, end))...
    ' and (5) context ' num2str(between_comparisons(5, 1)) ' and ' num2str(between_comparisons(5, end))...
    ' and (6) context ' num2str(between_comparisons(6, 1)) ' and ' num2str(between_comparisons(6, end))])
elseif size(between_comparisons,1) == 7
display(['between comparison is average of all distances between (1) context '...
    num2str(between_comparisons(1, 1)) ' and ' num2str(between_comparisons(1, end))...
    ' and (2) context ' num2str(between_comparisons(2, 1)) ' and ' num2str(between_comparisons(2, end))...
    ' and (3) context ' num2str(between_comparisons(3, 1)) ' and ' num2str(between_comparisons(3, end))...
    ' and (4) context ' num2str(between_comparisons(4, 1)) ' and ' num2str(between_comparisons(4, end))...
    ' and (5) context ' num2str(between_comparisons(5, 1)) ' and ' num2str(between_comparisons(5, end))...
    ' and (6) context ' num2str(between_comparisons(6, 1)) ' and ' num2str(between_comparisons(6, end))...
    ' and (7) context ' num2str(between_comparisons(7, 1)) ' and ' num2str(between_comparisons(7, end))])

end
    
    