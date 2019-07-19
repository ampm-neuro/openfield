function [all_within, all_between] = ALL_heatmaps_or_HD_firstwhite
%runs place cell function on every cell

count = 0;
pop = 0;

all_within = [];
all_between = [];

%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);

%ITERATE THROUGH SUBJECTS
%
for subject = 3%1:length_subjects
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
        for session = length_sessions
            day_file = file_list_sessions(session,1).name
            
            %load session
            load(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\' num2str(current_session_type) '\' num2str(day_file)]);

            if isempty(clusters)
                continue
            elseif size(clusters,1)<8
                %continue
            end
            
            pop > pop+1;
            if pop>1
                return
            end
           
            %INSERT FUNCTION HERE
            for cf = clusters(:,1)'

                %spatial heat maps
                %{
                bins = 50;
                figure; axis off; hold on
                [~,~,~,rm] = rate_mtx(eptrials(eptrials(:,6)==3,:), cf, bins); %skaggsmoothed rm
                rm = inpaint_nans(smooth2a(rm,1));
                imagesc(rm); 
                colormap jet
                colorbar
                title([num2str(current_rat) '\' num2str(current_session_type) '\' num2str(day_file) '\' num2str(cf)])
                set(gca,'TickLength',[0, 0]);
                axis square
                %}
                
                

            end
            
            %hd line plots
            %
            rate_dist = cell(4,1);
            
            %open new figure for each neurons
            for cf = 1:size(clusters,1)
               figure; 
            end
            
            if ~isequal(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))', 1:4) &&...
                    ~isequal(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))', 5:8)
                continue
            end

            for ctx = unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))'
                [rate_distributions] = hd_cell(eptrials(eptrials(:,6)==ctx,:), clusters(:,1), 11, 10, 1);              
                rate_dist{ctx} = rate_distributions;
            end
            
            %empty contexts
            if ismember(1, unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6)))
                withins = [1 2; 3 4];
                betweens = [1 3; 1 4; 2 3; 2 4];
                within_corrs = nan(size(clusters,1),size(withins,1));
                between_corrs = nan(size(clusters,1),size(betweens,1));
                for c = 1:size(clusters,1)
                    for wc = 1:size(withins,1)
                        within_corrs(c, wc) = corr(rate_dist{withins(wc,1)}(c,:)', rate_dist{withins(wc,2)}(c,:)');
                    end
                    for bc = 1:size(betweens,1)
                        between_corrs(c, bc) = corr(rate_dist{betweens(bc,1)}(c,:)', rate_dist{betweens(bc,2)}(c,:)');
                    end
                end
                
            %object contexts
            elseif ismember(5, unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6)))
                withins = [5 6];
                betweens = [7 8];
                within_corrs = nan(size(clusters,1),size(withins,1));
                between_corrs = nan(size(clusters,1),size(betweens,1));
                for c = 1:size(clusters,1)
                    for wc = 1:size(withins,1)
                        within_corrs(c, wc) = corr(rate_dist{withins(wc,1)}(c,:)', rate_dist{withins(wc,2)}(c,:)');
                    end
                    for bc = 1:size(betweens,1)
                        between_corrs(c, bc) = corr(rate_dist{betweens(bc,1)}(c,:)', rate_dist{betweens(bc,2)}(c,:)');
                    end
                end
            end
            
            %average across comparisons
            within_corrs = mean(within_corrs, 2);
            between_corrs = mean(between_corrs, 2);

            %load output
            all_within = [all_within; within_corrs];
            all_between = [all_between; between_corrs];

        end
    end
end


end