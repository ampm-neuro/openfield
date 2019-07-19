function ALL_generic
%runs place cell function on every cell

count = 0;

%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);

%ITERATE THROUGH SUBJECTS
%
for subject = 1:length_subjects
    current_rat = file_names_subjects{:}(subject,1).name
    file_list_session_type = dir(['neurodata/' num2str(current_rat) '/']);
    file_list_session_type(1:2) = [];
    file_names_session_types = {file_list_session_type([file_list_session_type(:).isdir])};
    length_session_types = size(file_names_session_types{:}, 1);
    
    %ITERATE THROUGH SESSION TYPE
    %
    for session_type = 1:length_session_types
        current_session_type = file_names_session_types{:}(session_type,1).name;
        file_list_sessions = dir(['neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/*.mat']);
        length_sessions = size(file_list_sessions,1);
        
        %ITERATE THROUGH SESSIONS
        %
        for session = 1:length_sessions
            day_file = file_list_sessions(session,1).name;
            
            %load session
            load(['/Users/ampm/Documents/MATLAB/PCia/neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/' num2str(day_file)]);

            %INSERT FUNCTION HERE
            
            
            if isempty(clusters)
                continue
            end
           
            
            for clust = unique(clusters(:,1))'
                
                
                plot_FR_by_col(eptrials, clust, 7, 30, 500, 500, 0)
                
                
            end
            
            
        
            
        end
    end
end


end