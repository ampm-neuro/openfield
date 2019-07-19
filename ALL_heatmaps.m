function [rms, hd_props, rat_num] = ALL_heatmaps
%runs place cell function on every cell

clust_count = 0;
bins = 40;
rms = nan(4,bins^2);
context_titles{1} = 'B1';
context_titles{2} = 'B2';
context_titles{3} = 'W1';
context_titles{4} = 'W2';
context_titles{5} = 'A1';
context_titles{6} = 'A2';
context_titles{7} = 'O1';
context_titles{8} = 'O2';
rat_num = [];


%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);
hd_props = nan(4, bins, bins, 360, 161);

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
    for session_type = 1%:length_session_types
        current_session_type = file_names_session_types{:}(session_type,1).name
        file_list_sessions = dir(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\' num2str(current_session_type) '\*.mat']);
        %any2csv(file_list_sessions, '|', 1)
        length_sessions = size(file_list_sessions,1);
        
        %ITERATE THROUGH SESSIONS
        %
        for session = 1:length_sessions
            day_file = file_list_sessions(session,1).name
            
            %load session
            load(['C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\' num2str(current_rat) '\' num2str(current_session_type) '\' num2str(day_file)]);
            
            %{
            if ~isequal(unique(eptrials(~isnan(eptrials(:,6)) & eptrials(:,6)<9,6))', 5:8)
                continue
            end
            %}
            
            if isempty(clusters)
                continue
            end
            %clusters = clusters(clusters(:,2)>2,:);
            %if isempty(clusters)
            %    continue
            %end
            
           
            %INSERT FUNCTION HERE
            for cf = clusters(:,1)'
                %figure; axis off; hold on
                %ha = tight_subplot(2,2,[.05 -.3],[.005 .005],[.005 .005]);
                clust_count = clust_count+1;
                
                vis_count = 0;
                for visit = 1:4
                    vis_count = vis_count+1;
                    %current visit
                    if  mean(eptrials(~isnan(eptrials(:,6)) &...
                        eptrials(:,6)<9,6)) < 4.5
                        visit_idx = eptrials(:,6)==visit; %BBWW
                    else
                        visit_idx = eptrials(:,6)==visit+4;
                    end
                    
                    %plot heatmap
                    %subplot(1,4, vis_count);
                    %axes(ha(vis_count))
                    [rm, ~, ~, ~, HD_proportions] = rate_mtx(eptrials(visit_idx,:), cf, bins);
                    %rm = smooth2a(rm,1);
                    rms(visit,:,clust_count) = rm(:);
                    hd_props(visit, :,:,:, clust_count) = HD_proportions;
                    rat_num = [rat_num; subject];
                end
            end
        end
    end
end


end