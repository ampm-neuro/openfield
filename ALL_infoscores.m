function iscores = ALL_infoscores(bw_pop_cell_idx)
%calculates information score for every cell

count = 0;

iscores = [];

%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);

%ITERATE THROUGH SUBJECTS
%
for subject = 1:length_subjects
    current_rat = file_names_subjects{:}(subject,1).name;
    file_list_session_type = dir(['neurodata/' num2str(current_rat) '/']);
    file_list_session_type(1:2) = [];
    file_names_session_types = {file_list_session_type([file_list_session_type(:).isdir])};
    length_session_types = size(file_names_session_types{:}, 1);
    
    %ITERATE THROUGH SESSION TYPE
    %
    for session_type = 1
        current_session_type = file_names_session_types{:}(session_type,1).name;
        file_list_sessions = dir(['neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/*.mat']);
        length_sessions = size(file_list_sessions,1);
        
        %ITERATE THROUGH SESSIONS
        %
        for session = 1:length_sessions
            day_file = file_list_sessions(session,1).name;
            
            %load session
            load(['neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/' num2str(day_file)]);

            %INSERT FUNCTION HERE
            if isempty(clusters)
                continue
            end
            
            %constrain clusters
            cluster_confidence = [3 4 5];
            cluster_region = [0 1 2];
            cluster_idx = ismember(clusters(:,2), cluster_confidence) & ismember(clusters(:,4), cluster_region);
            clusters = clusters(cluster_idx, :);
                
            %skip sessions without cells
            if isempty(clusters)
                continue
            end
           
            %context of interest
            coi = 3;
            
            %calculate scores
            %if mean(bw_pop_cell_idx(count+1:count+length(clusters(:,1)))) == 5
                is_out = info_score(eptrials(eptrials(:,6)==coi,:), 20, clusters(:,1)');
                iscores = [iscores; is_out'];
            %end
            
            count = count+length(clusters(:,1));
            
            
        
            
        end
    end
end


end