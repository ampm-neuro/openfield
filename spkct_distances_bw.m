function [within_comps, between_comps] = spkct_distances_bw(spike_counts, contexts, within_comparisons, between_comparisons, method1, method2, time_window)
% distances(spike_counts, contexts, comparisons, method1, method2)
%
% calculates distances between HD clouds of spikecounts computed with 
% ALL_countwindow
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
        within_comps = nan(size(within_comparisons,1), 1);
        for wc = 1: size(within_comps,1)
            within_comps(wc) = pdist([ctxt_means(within_comparisons(wc,1), :); ctxt_means(within_comparisons(wc,2), :)], 'euclidean')./dc;
        end

        %betweens
        between_comps = nan(size(between_comparisons,1), 1);
        for bc = 1: size(between_comps,1)
            between_comps(bc) = pdist([ctxt_means(between_comparisons(bc,1), :); ctxt_means(between_comparisons(bc,2), :)], 'euclidean')./dc;
        end
        
        %figure
        figure; 
        bar([1 2], [mean(between_comps) mean(within_comps)])
        xlim([.5 2.5])
        xticks(1:2)
        xticklabels({'Between', 'Within'})
        xlabel contexts
        ylabel('euclid distance')
        
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
        within_comps = nan(size(spike_counts,1),1);
        between_comps = nan(size(spike_counts,1),1);
        
        %all samples
        allsamps = 1:size(spike_counts,1);
        
        %for each sample
        for samp = allsamps
            
            %spike_counts sans samp
            spike_counts_dropone = spike_counts(allsamps~=samp,:);
            contexts_dropone = contexts(allsamps~=samp);
            
            %recalculate means for each samp (dropping test samp)
            ctxt_means_dropone = nan(length(rlvt_ctxts), size(spike_counts_dropone,2));
            for cm = 1:size(ctxt_means_dropone,1)
                ctxt_means_dropone(cm, :) = mean(spike_counts_dropone(contexts_dropone == rlvt_ctxts(cm),:));
            end

            %calculate dist from samp to every mean (also prod. extraneous)
            dist_to_means = pdist([spike_counts(samp,:); ctxt_means_dropone]); dist_to_means = dist_to_means(1:size(ctxt_means_dropone,1))./dc;
            
            %dist to self
            dts = dist_to_means(rlvt_ctxts == contexts(samp));
            
            %current withins and betweens
            cc_w = within_comparisons(logical(sum(within_comparisons == contexts(samp),2)),:); 
            cc_w = cc_w(cc_w~=contexts(samp));
            %cc_w = mean(cc_w);
            cc_b = between_comparisons(logical(sum(between_comparisons == contexts(samp),2)),:); 
            cc_b = cc_b(cc_b~=contexts(samp));

            %assign the within distance
            within_comps(samp) = mean(dist_to_means(cc_w)) -1;

            %assign the between distances
            between_comps(samp) = mean(dist_to_means(cc_b)) -1 ;
            
        end
        
        %remove nans
        within_comps = within_comps(~isnan(within_comps));
        between_comps = between_comps(~isnan(between_comps));
    
        %figure
        %

        %average distance from samp to contexts of the same type
        cell_e{1} = mean(between_comps,2);% - ones(size(between_comps(:,1)));

        %average distance from samp to contexts of the opposite type
        cell_e{2} = nanmean(within_comps,2);% - ones(size(within_comps(:,1)));

        %the above means combine 4 categories with distinct distances (not
        %perfectly related to column. plot lines of columns to see.

        %plot bar with errorbars
        errorbar_barplot(cell_e);
        xticks(1:2)
        xticklabels({'Between', 'Within'})
        xlabel contexts
        ylabel('distance')
    
    
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
        within_comps = nan(size(within_comparisons,1), 1);
        for wc = 1: size(within_comps,1)
            %average forward and reverse
            mahal_forward_w = mahal(ctxt_means(within_comparisons(wc,1), :), ctxt_clouds{within_comparisons(wc,2)})./dc;
            mahal_backward_w = mahal(ctxt_means(within_comparisons(wc,2), :), ctxt_clouds{within_comparisons(wc,1)})./dc;
            within_comps(wc) = mean([mahal_forward_w mahal_backward_w]);
        end

        %betweens
        between_comps = nan(size(between_comparisons,1), 1);
        for bc = 1: size(between_comps,1)            
            %average forward and reverse
            mahal_forward_b = mahal(ctxt_means(between_comparisons(bc,1), :), ctxt_clouds{between_comparisons(bc,2)})./dc;
            mahal_backward_b = mahal(ctxt_means(between_comparisons(bc,2), :), ctxt_clouds{between_comparisons(bc,1)})./dc;
            between_comps(bc) = mean([mahal_forward_b mahal_backward_b]);
        end
        
        %figure
        figure; 
        bar([1 2], [mean(between_comps) mean(within_comps)])
        xlim([.5 2.5])
        xticks(1:2)
        xticklabels({'Between', 'Within'})
        xlabel contexts
        ylabel('mahal distance')
        
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
        within_comps = nan(size(spike_counts,1),2);
        between_comps = nan(size(spike_counts,1),2);
        
        %all samples
        allsamps = 1:size(spike_counts,1);
        
        %for each sample
        for samp = allsamps
            
            %current sample spikecounts
            c_samp = spike_counts(samp,:);
            
            %spike_counts sans samp
            spike_counts_dropone = spike_counts(allsamps~=samp,:);
            contexts_dropone = contexts(allsamps~=samp);
            
            %preallocate
            cloud_dists = nan(length(rlvt_ctxts),1);
            
            %for each other context
            for cds = 1:length(cloud_dists)
                %{
                %current context cloud
                c_ctxt = spike_counts_dropone(contexts_dropone == rlvt_ctxts(cds),:);
                
                %distance from samp to cloud
                cloud_dists(cds) = mahal(c_samp, c_ctxt)./dc;
                %}
                
                %distance from samp to cloud
                cloud_dists(cds) = mahal(spike_counts(samp,:), spike_counts_dropone(contexts_dropone == rlvt_ctxts(cds),:))./dc;


            end
            
            %current withins and betweens
            cc_w = within_comparisons(logical(sum(within_comparisons == contexts(samp),2)),:); 
            %cc_w = cc_w(cc_w~=contexts(samp));
            cc_b = between_comparisons(logical(sum(between_comparisons == contexts(samp),2)),:); 
            cc_b = cc_b(cc_b~=contexts(samp));
            
            %assign the within distance
            within_comps(samp,:) = mean(cloud_dists(cc_w));%mean(cloud_dists(cc_w));

            %assign the between distances
            between_comps(samp,:) = cloud_dists(cc_b);%mean(cloud_dists(cc_b));
            
        end
        
        
        %figure
        %
        
        %average distance from samp to contexts of the same type
        cell_e{1} = mean(between_comps,2);
        
        %average distance from samp to contexts of the opposite type
        cell_e{2} = nanmean(within_comps,2);
        
        %the above means combine 4 categories with distinct distances (not
        %perfectly related to column. plot lines of columns to see.
        
        %plot bar with errorbars
        errorbar_barplot(cell_e);
        xticks(1:2)
        xticklabels({'Between', 'Within'})
        xlabel contexts
        ylabel('mahal distance')

        %remove nans
        %within_comps = within_comps(~isnan(within_comps));
        %between_comps = between_comps(~isnan(between_comps));
        
        
    end
end

end
    
    