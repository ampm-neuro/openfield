function [rate_dists, corrs] = ALL_HD_ctxt(varargin)
%runs hd_rel_obj function on every cell, context, and ~quarter map

if nargin>0
    rate_dists = varargin;
    corrs = varargin;
end

pop = 0;

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
    for session_type = 1
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
            
            if ~isequal(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))', 1:4)
                continue
            end

            if exist(['C:\Users\ampm1\Desktop\context paper\HD_bw\'...
                    'HD_bw_' num2str(current_rat) '_' num2str(session) '_' ...
                    num2str(floor(clusters(end,1))) '_' ...
                    num2str((clusters(end,1)-floor(clusters(end,1)))*10) '.pdf'], 'file')
                continue
            end

            
           
            %INSERT FUNCTION HERE
            %open new figure for each neurons
            for cf = 1:size(clusters,1)
               figure; hold on;
                %preallocateish
                if exist('rate_dists','var')
                    rate_dists{length(rate_dists) + 1} = cell(1,4);
                else
                    rate_dists{cf} = cell(1,4);
                    corrs = [];
                end
            end
            num_figs = size(findobj('Type', 'figure'),1);
            target_figs = (num_figs-size(clusters,1)+1):num_figs;

            for ctx = unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))'

                %calculate zscored rate distribution
                rd_all = hd_cell(eptrials(eptrials(:,6)==ctx,:), clusters(:,1), 11, 7, 1);   
                rd_all = rd_all-repmat(nanmean(rd_all, 2), 1, size(rd_all, 2));
                rd_all = rd_all./repmat(nanstd(rd_all, [], 2), 1, size(rd_all,2));
                
                %plot on each cell's figure
                for coi = 1:size(clusters,1)
                    figure(target_figs(coi)); hold on
                    
                    if ismember(ctx, 1:2)
                        plot(rd_all(coi,:), 'k', 'linewidth', 1.1, 'color', [0    0.5    0.7]) %blue for black contexts
                    elseif ismember(ctx, 3:4)
                        plot(rd_all(coi,:), 'k', 'linewidth', 1.1, 'color', [0.9    0.3    0.1]) %orange for white contexts
                    
                        if ctx==4
                            set(gca,'TickLength',[0, 0]);
                            box off; xlim([0 360]); ylim([-5 5])
                            xticks([1 90 180 270 360])
                            xticklabels({'180', '270', '0', '90', '179'})
                            xlabel('Head direction')
                            ylabel('FR (z Hz)')
                        end
                    end

                    %load distribution output
                    rate_dists{length(rate_dists) - size(clusters,1) + coi}{1,ctx} = rd_all(coi,:);
                    
                end
            end
            
            
            %compute correlations between context distributions
            if ~isnan(mean(mean([rate_dists{coi}{1,1}; rate_dists{coi}{1,2};...
                    rate_dists{coi}{1,3}; rate_dists{coi}{1,4}])))
                
                within_corr = mean([...
                        corr(rate_dists{coi}{1,1}', rate_dists{coi}{1,2}') ...
                        corr(rate_dists{coi}{1,3}', rate_dists{coi}{1,4}')]);
                between_corr = mean([...
                        corr(rate_dists{coi}{1,1}', rate_dists{coi}{1,3}') ...
                        corr(rate_dists{coi}{1,1}', rate_dists{coi}{1,4}') ...
                        corr(rate_dists{coi}{1,2}', rate_dists{coi}{1,3}') ...
                        corr(rate_dists{coi}{1,2}', rate_dists{coi}{1,4}')]);
                corrs = [corrs; [within_corr between_corr]];
            
            else
                corrs = [corrs; [nan nan]];
            end
            
            
            %save each figure
            %
            for ic = 1:length(target_figs)
               figure(target_figs(ic))
               var_name = ['HD_bw' num2str(current_rat) '_' num2str(session) '_' ...
                   num2str(floor(clusters(ic,1))) '_' num2str((clusters(ic,1)-floor(clusters(ic,1)))*10)];%, ...
                   %'_hem' num2str() '_rg' num2str()]; 
               print(['C:\Users\ampm1\Desktop\context paper\HD_bw\' var_name],...
                   '-dpdf', '-painters', '-bestfit')
               close
            end
            close all
            %}
            save('HD_bw_distribs')
        end
    end
end

end