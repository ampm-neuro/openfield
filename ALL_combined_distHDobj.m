function [all_rd_heatmaps] = ALL_combined_distHDobj(varargin)
%runs hd_rel_obj function on every cell, context, and ~quarter map

%edges
HD_edges = 0:10:359.000000001;
dist_edges = 0:0.025:1;

count_c=0;

if nargin>0
    rate_dists = varargin;
end

%preallocate
all_rd_heatmaps = [];

%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);

%ITERATE THROUGH SUBJECTS
%
for subject = 1:length_subjects
    current_rat = file_names_subjects{:}(subject,1).name
    file_list_session_type = dir(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\']);
    file_list_session_type(1:2) = [];
    file_names_session_types = {file_list_session_type([file_list_session_type(:).isdir])};
    length_session_types = size(file_names_session_types{:}, 1);
    
    %ITERATE THROUGH SESSION TYPE
    %
    for session_type = 2
        current_session_type = file_names_session_types{:}(session_type,1).name;
        file_list_sessions = dir(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\' num2str(current_session_type) '\*.mat']);
        %any2csv(file_list_sessions, '|', 1)
        length_sessions = size(file_list_sessions,1);
        
        %ITERATE THROUGH SESSIONS
        %
        for session = 1:length_sessions
            day_file = file_list_sessions(session,1).name
            
            %load session
            load(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\' num2str(current_session_type) '\' num2str(day_file)]);

            if isempty(clusters)
                continue
            elseif size(clusters,1)<8
                %continue
            end
            
            if ~isequal(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))', 5:8)
                continue
            end
            
            %INSERT FUNCTION HERE
            count_c_hold = count_c;
            for ctx = unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))'...
                    - (min(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6)))-1)
                objpos = object_grid(ctx+4, 1, 0);
                

                %allocentric
                occupancy_hd = eptrials(eptrials(:,4)==1, 9);
                occupancy_hd = histcounts(occupancy_hd, HD_edges)./100;
                for coi = 1:size(clusters,1)
                    count_c = count_c+1;
                    spike_counts_hd = eptrials(eptrials(:,4)==clusters(coi,1) & eptrials(:,6)==ctx+4, 9);
                    all_rd_heat(count_c, :, ctx) = histcounts(spike_counts_hd, HD_edges)./occupancy_hd;
                end
                if ismember(ctx,1:3)
                    count_c = count_c_hold;
                else
                    count_c = count_c_hold+size(clusters,1);
                end
                
                for obj = 1:4
                    objpos_local = objpos(obj,:);
                    
                    %calculate hd relative to obj while the rat is occupying
                    %each of 4 spatial areas around the object
                    rm_NW = objfunctions(eptrials, ctx+4, clusters, obj, [0 objpos_local(1) objpos_local(2) 1], HD_edges, dist_edges);
                        rm_NW(isinf(rm_NW))=nan;
                    rm_NE = objfunctions(eptrials, ctx+4, clusters, obj, [objpos_local(1) 1 objpos_local(2) 1], HD_edges, dist_edges);
                        rm_NE(isinf(rm_NE))=nan;
                    rm_SE = objfunctions(eptrials, ctx+4, clusters, obj, [objpos_local(1) 1  0 objpos_local(2)], HD_edges, dist_edges);
                        rm_SE(isinf(rm_SE))=nan;
                    rm_SW = objfunctions(eptrials, ctx+4, clusters, obj, [0 objpos_local(1) 0 objpos_local(2)], HD_edges, dist_edges);
                        rm_SW(isinf(rm_SW))=nan;
                    
                    for coi = 1:size(clusters,1)
                        
                        %{
                        if coi == 5
                        rm_heatmap(rm_NW, coi); title(['Context ' num2str(ctx) ', object ' num2str(obj) ' NW'])
                        rm_heatmap(rm_NE, coi); title(['Context ' num2str(ctx) ', object ' num2str(obj) ' NE']) 
                        rm_heatmap(rm_SE, coi); title(['Context ' num2str(ctx) ', object ' num2str(obj) ' SE'])
                        rm_heatmap(rm_SW, coi); title(['Context ' num2str(ctx) ', object ' num2str(obj) ' SW'])
                        end
                        %}
                        
                        %{
                        if coi == 1
                            rm_heatmap(rm_NW, coi); title(['Context ' num2str(ctx) ', object ' num2str(obj) ' NW']);
                        end
                        %}
                        
                        mean_obj_rm(:, :, obj, ctx, coi) = nanmean(cat(3, rm_NW(:,:,coi), rm_NE(:,:,coi), rm_SE(:,:,coi), rm_SW(:,:,coi)), 3);
                    end
                    mean_obj_rm(isinf(mean_obj_rm))=nan;
                end
                
                %average across objects in the context
                for coi = 1:size(clusters,1)
                    mean_ctx_rm(:, :, ctx, coi) = nanmean(mean_obj_rm(:, :, :, ctx, coi),3);
                    mean_ctx_rm(isinf(mean_ctx_rm)) = nan;
                end
            end
            
            %average across contexts
            for coi = 1:size(clusters,1)
                mean_rm(:, :, :, coi) = nanmean(mean_ctx_rm(:, :, :, coi),3);
                mean_rm(isinf(mean_rm)) = nan;
                
                all_rd_heatmaps = cat(3, all_rd_heatmaps, mean_rm(:, :, coi)');
                
                %{
                rm_heatmap(mean_rm, coi)
                title(num2str(clusters(coi,1)))
                %}
                
                %save each figure
                %{
                var_name = ['DistHD_rel_obj_' num2str(current_rat) '_' num2str(session) '_' ...
                       num2str(floor(clusters(coi,1))) '_' num2str((clusters(coi,1)-floor(clusters(coi,1)))*10)];%, ...
                print(['C:\Users\ampm1\Desktop\context paper\DistHD_ratemaps\' var_name],...
                   '-dpdf', '-painters', '-bestfit')
                %}
            end
           
            
        end
    end
end

end

%calculates spike and occupancy counts (and ratematrix) from external functions
function rate_matrix = objfunctions(eptrials, ctx, clusters, obj, ba, HD_edges, dist_edges)

    %constrain eptrials
    eptrials = eptrials(eptrials(:,2)>ba(1) & eptrials(:,2)<ba(2) & ...
        eptrials(:,3)>ba(3) & eptrials(:,3)<ba(4) & eptrials(:,6)==ctx, :);
    
    %calculate counts for every cell
    occupancy_hd = eptrials(eptrials(:,4)==1, 9+obj);
    occupancy_dist = eptrials(eptrials(:,4)==1, 13+obj);
    occupancy_mtx = histcounts2(occupancy_dist, occupancy_hd, HD_edges, dist_edges)./100;
        
    %compute rate matrices for each cell   
    for coi = 1:size(clusters(:,1))
        spike_counts_hd = eptrials(eptrials(:,4)==clusters(coi,1), 9+obj);
        spike_counts_dist = eptrials(eptrials(:,4)==clusters(coi,1), 13+obj);
        spike_counts_mtx = histcounts2(spike_counts_dist, spike_counts_hd, HD_edges, dist_edges);
        
        %figure; imagesc(spike_counts_mtx); title spikes
        %figure; imagesc(occupancy_mtx); title pos
        
        
        rate_matrix(:, :, coi) = (spike_counts_mtx./occupancy_mtx);
    end
end

function rm_heatmap(rm, coi)
    rm_fig = rm(:, :, coi)';
    rm_fig = rm_fig(:,[floor(size(rm_fig,2)/2)+1:size(rm_fig,2) 1:floor(size(rm_fig,2)/2)]);
    rm_fig = smooth2a(rm_fig,1);
    figure; imagesc(rm_fig); 
    axis square; set(gca,'TickLength',[0, 0]); box off
    ylabel('Distance from object (m)')
    xlabel('Head direction (degrees from object)')
    xticks([1 size(rm_fig,2)/2 size(rm_fig,2)])
    xticklabels({'-180' 'Object' '180'});
    yticks([1 size(rm_fig,1)])
    yticklabels({'0' '1'});
    colorbar
end