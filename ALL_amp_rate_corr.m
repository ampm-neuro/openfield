function [corr_vals, FRs_amps, all_cell_idx, all_cell_idx_bw, all_cell_idx_ia, tt_id, all_amplitudes] = ALL_amp_rate_corr(subjects, sessions)
% iterate through folders on external hd and computer corr between
% amplitude and firing rate
%
% assign unique tt-session number to each amplitude r val

all_cell_idx = [];
all_cell_idx_bw = [];
all_cell_idx_ia = [];
corr_vals = [];
FRs_amps = [];
tt_id = []; %session number + (tt number / 100)
all_amplitudes = [];

idx_low = 0;
idx_hi = 0;

%path to raw neuralynx files (origin)
%path_origin = '/Volumes/LaCie/PCAlte/';
%path_origin = '/Volumes/ampm_PC2Mac/';
path_origin = 'D:\array_5_15_2017\Project\PC Identity Arrangement\';

%get all the things in neurodata folder...
file_list_subjects = dir(path_origin);

%hard coded erasure of irrelevant directory folders
file_list_subjects = file_list_subjects(4:7);

%exclude non-folders
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%number of folders
length_subjects = size(file_names_subjects{:},1);


%iterate through subjects
all_sesh_num = 0;
for subject = subjects%1:length_subjects

    %print update
    rat = file_names_subjects{:}(subject,1).name

    %get all the things in subject folder...
    file_list_sessions = dir(strcat(path_origin, num2str(file_names_subjects{:}(subject,1).name)));

    %hard coded erasure from file_list_stages of irrelevant directory folders
    file_list_sessions(1:2) = [];

    %exclude non-folders
    file_names_sessions = {file_list_sessions([file_list_sessions(:).isdir])};
    
    %number of folders
    length_sessions = size(file_names_sessions{:},1);
    
    %iterate through session folders
    for session = sessions%1:length_sessions
        all_sesh_num = all_sesh_num + 1;
        
        %current session file name
        session_file = file_names_sessions{:}(session,1).name
        session_number = num2str(session);
        
        if isempty(strfind(session_file, 'rec'))
            continue
        end
        
        %identify stage from session file name
        if ~isempty(strfind(session_file, 'a1'))
            stage = 'object_arrangement';
        elseif ~isempty(strfind(session_file, 'A1'))
            stage = 'object_arrangement';
        else
            stage = 'black_and_white';
        end 
            

       
        %load session
        filename = strcat(path_origin, rat, '\', session_file);
        
        %load excel files
        clusters = xlsread([filename '\Cluster_descriptions.xlsx']);
        if isequal(size(clusters), [2,1])
            clusters = [];
            continue
        else
            clusters = clusters(~isnan(clusters(:,1)),:);
        end
        
        %constrain clusters
        cluster_confidence = [3 4 5];
        cluster_region = [0 1 2];
        cluster_idx = ismember(clusters(:,3), cluster_confidence) & ismember(clusters(:,5), cluster_region);
        clusters = clusters(cluster_idx, :);
        if isempty(clusters)
            continue
        end
                  
        %GET TT        
        [cv_local, fa_local, ~, ~, all_cell_idx_local, amplitudes] = GetTT_waveforms(filename, cluster_idx, session_file);    
        
        %FIX
        %cluster_idx = cluster_idx(ismember(clusters(:,3), cluster_confidence));

        %load
        if sum(cluster_idx) > 0
            
            all_cell_idx = [all_cell_idx; all_cell_idx_local];
            if strcmp(stage, 'black_and_white')
                all_cell_idx_bw = [all_cell_idx_bw; all_cell_idx_local];
            elseif strcmp(stage, 'object_arrangement')
                all_cell_idx_ia = [all_cell_idx_ia; all_cell_idx_local];
            end
            corr_vals = [corr_vals; cv_local];
            FRs_amps = [FRs_amps; fa_local];
            
            tt_id = [tt_id; clusters]; 
            all_amplitudes = [all_amplitudes; amplitudes];
            
        end
        
    end

end

end