
load('all_amp_drift')

if ~exist('all_amplitudes','var')
    corr_vals=[];
    FRs_amps=[];
    all_cell_idx=[]; 
    all_cell_idx_bw=[]; 
    all_cell_idx_ia=[]; 
    tt_id=[]; 
    all_amplitudes=[];
    subj_idx = [];
    sesh_idx = [];
    count_hold = 0;
end

count = 1;

for subject = 1:4 
    

    %identify number of sessions from folders
    path_origin = 'D:\array_5_15_2017\Project\PC Identity Arrangement\';
    file_list_subjects = dir(path_origin);
    file_list_subjects = file_list_subjects(4:7);
    file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
    rat = file_names_subjects{:}(subject,1).name;

    %get all the things in subject folder...
    file_list_sessions = dir(strcat(path_origin, num2str(file_names_subjects{:}(subject,1).name)));
    file_list_sessions(1:2) = [];
    file_names_sessions = {file_list_sessions([file_list_sessions(:).isdir])};
    length_sessions = size(file_names_sessions{:},1);
    
    for sesh = 1:length_sessions
    
        if count <=count_hold
            count = count+1;
            continue
        end
        
        %current session file name
        session_file = file_names_sessions{:}(sesh,1).name
        session_number = num2str(sesh);
        
        if ~contains(session_file, 'rec')
            continue
        end

        %identify stage from session file name
        if contains(session_file, 'a1')
            stage = 'object_arrangement';
        elseif contains(session_file, 'A1')
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

        %load
        if sum(cluster_idx) > 0

            if strcmp(stage, 'black_and_white')
                all_cell_idx_bw_sesh = all_cell_idx_local;
                all_cell_idx_bw=[all_cell_idx_bw; all_cell_idx_bw_sesh];
            elseif strcmp(stage, 'object_arrangement')
                all_cell_idx_ia_sesh = all_cell_idx_local;
                all_cell_idx_ia=[all_cell_idx_ia; all_cell_idx_ia_sesh];
            end
        
            corr_vals=[corr_vals; cv_local];
            FRs_amps=[FRs_amps; fa_local];
            all_cell_idx=[all_cell_idx; all_cell_idx_local]; 
            tt_id=[tt_id; clusters]; 
            all_amplitudes=[all_amplitudes; amplitudes];

            subj_idx = [subj_idx; repmat(subject, size(amplitudes))];
            sesh_idx = [sesh_idx; repmat(sesh, size(amplitudes))];
        
        end
        
        
        count = count+1;
        count_hold = count_hold+1

        save('all_amp_drift')
    
    end
    
end