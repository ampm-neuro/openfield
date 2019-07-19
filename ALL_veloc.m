function [veloc_mtx_visitorder, veloc_mtx_sortctxts] = ALL_veloc
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

veloc_mtx_visitorder = [];
veloc_mtx_sortctxts = [];
%cumdist_mtx = [];
%veloc_mtx_first = [];
%cumdist_mtx_first = [];
%veloc_mtx_second = [];
%cumdist_mtx_second = [];

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
    for session_type = 1%:length_session_types
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
            [velocs] = ctxt_vels(eptrials, context_ids);

            veloc_mtx_visitorder = [veloc_mtx_visitorder; velocs];

            [~, srt_idx] = sort(context_ids);
            veloc_mtx_sortctxts = [veloc_mtx_sortctxts; velocs(srt_idx)];
            
        
            
        end
    end
end






end

