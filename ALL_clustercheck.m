function ALL_clustercheck
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
            
            
            if isempty(clusters) && isempty(unique(eptrials(~isnan(eptrials(:,4)) & eptrials(:,4)>1, 4)))
                continue
            end
           
            %clust lists
            c_excel = clusters(:,1);
            c_eptrials = unique(eptrials(~isnan(eptrials(:,4)) & eptrials(:,4)>1, 4));
            
            %if inequal
            if ~isequal(c_excel, c_eptrials)
                
                %show inequal file and clust lists
                [num2str(current_rat) '/' num2str(current_session_type) '/' num2str(day_file)]
                c_excel = clusters(:,1)
                c_eptrials = unique(eptrials(~isnan(eptrials(:,4)) & eptrials(:,4)>1, 4))
                
                %if same length
                if length(c_excel) == length(c_eptrials)
                    
                    %find inequal elements
                    bad_elem = [c_excel(c_excel~=c_eptrials) c_eptrials(c_excel~=c_eptrials)]
                    
                %if different lengths
                else
                    
                    %find missing elements
                    if length(c_excel) > length(c_eptrials)
                        
                        missing_elem = c_excel(~ismember(c_excel, c_eptrials));
                        bad_elem = [missing_elem nan(size(missing_elem))]
                        
                    else
                        
                        missing_elem = c_eptrials(~ismember(c_eptrials, c_excel));
                        bad_elem = [nan(size(missing_elem)) missing_elem]
                        
                    end
                    
                end
                
                
                
            end
            
        
            
        end
    end
end


end