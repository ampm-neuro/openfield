function rate_dists = ALL_dist_relative_to_obj(varargin)
%runs hd_rel_obj function on every cell, context, and ~quarter map

if nargin>0
    rate_dists = varargin;
end
colors = get(gca,'ColorOrder');close

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

            %{
            if exist(['C:\Users\ampm1\Desktop\context paper\distS_rel_obj\'...
                    'HD_rel_obj_' num2str(current_rat) '_' num2str(session) '_' ...
                    num2str(floor(clusters(end,1))) '_' ...
                    num2str((clusters(end,1)-floor(clusters(end,1)))*10) '.pdf'], 'file')
                continue
            end
            %}

            
           
            %INSERT FUNCTION HERE
            %open new figure for each neurons
            for cf = 1:size(clusters,1)
               figure; set(gcf, 'position', [500 200 1000 800])
               
                %preallocateish
                if exist('rate_dists','var')
                    rate_dists{length(rate_dists) + 1} = cell(5,4);
                    spike_counts{length(spike_counts) + 1} = cell(5,4);
                    occupancy_counts{length(occupancy_counts) + 1} = cell(5,4);
                    rate_distr_unsm{length(rate_distr_unsm) + 1} = cell(5,4);
                else
                    rate_dists{cf} = cell(5,4);
                    spike_counts{cf} = cell(5,4);
                    occupancy_counts{cf} = cell(5,4);
                    rate_distr_unsm{cf} = cell(5,4);
                end
            end
            num_figs = size(findobj('Type', 'figure'),1);
            target_figs = (num_figs-size(clusters,1)+1):num_figs;

            for ctx = unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))'
                objpos = object_grid(ctx, 1, 0);

                for obj = 1:4
                    objpos_local = objpos(obj,:);
                    
                    %calculate dist relative to obj 
                    [rate_distributions, spike_count_distr, dist_counts, rate_distr_unsmoothed] = dist2obj_rate(eptrials(eptrials(:,6)==ctx,:), clusters(:,1), 11, 3, 1, 9+obj, [0 objpos_local(1) objpos_local(2) 1]);
                                        
                    for coi = 1:size(rate_distributions,1)
                        %load output
                        rate_dists{length(rate_dists) - size(clusters,1) + coi}{obj,ctx-4} = rate_distributions(coi,:);
                        spike_counts{length(spike_counts) - size(clusters,1) + coi}{obj,ctx-4} = spike_count_distr(coi,:);
                        occupancy_counts{length(occupancy_counts) - size(clusters,1) + coi}{obj,ctx-4} = dist_counts';
                        rate_distr_unsm{length(rate_distr_unsm) - size(clusters,1) + coi}{obj,ctx-4} = rate_distr_unsmoothed(coi,:);
                        
                        %load cell's figure
                        figure(target_figs(coi))
                        
                        %add object rows
                        if obj == 1
                            subplot(6,4, ctx-4)
                            trlfree_heatmap(eptrials(eptrials(:,6)==ctx,:), clusters(coi,1), 40, 1, ctx);
                            %object_grid(ctx);
                            yticks([])
                            xticks([])
                            box on
                        end
                        
                        %zscore rate distribution
                        %{
                        rd = rate_distributions(coi,:)-repmat(nanmean(rate_distributions(coi,:)),size(rate_distributions(coi,:)));
                        rd = rd/nanstd(rd);

                        %plot
                        subplot(6,4, 4*(obj-1)+(ctx-4) + 4)
                        plot(rd, 'color', colors(obj,:), 'linewidth', 1.1)
                        set(gca,'TickLength',[0, 0]); 
                        box off; xlim([0 size(rd,2)]); ylim([-5 5])
                        if ctx~=5
                            yticks([])
                        else
                            ylabel('FR (z Hz)')
                        end
                        
                        
                        subplot(6,4, 20+(ctx-4))
                        hold on
                        plot(rd, 'color', colors(obj,:), 'linewidth', 1.1)
                        set(gca,'TickLength',[0, 0]); 
                        box off; xlim([0 size(rd,2)]); ylim([-5 5])
                        
                        
                        %xticks([1 size(rd,2)/2 size(rd,2)])
                        %xticklabels({'-180', ['O' num2str(obj)], '180'})
                        if ctx~=5
                            yticks([])
                        else
                            ylabel('FR (z Hz)')
                        end
                        %}
                    end
                end
            end
            
            %save each figure
            %{
            for ic = 1:length(target_figs)
               figure(target_figs(ic))
               var_name = ['dist_rel_obj_' num2str(current_rat) '_' num2str(session) '_' ...
                   num2str(floor(clusters(ic,1))) '_' num2str((clusters(ic,1)-floor(clusters(ic,1)))*10)];%, ...
                   %'_hem' num2str() '_rg' num2str()]; 
               print(['C:\Users\ampm1\Desktop\context paper\dist_rel_obj\' var_name],...
                   '-dpdf', '-painters', '-bestfit')
               close
            end
            %}
            close all
            save('HD_dist_distribs_counts')
        end
    end
end

end