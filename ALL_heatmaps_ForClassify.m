function [out_mtx, out_mtx_z] = ALL_heatmaps_ForClassify(bins)
%runs place cell function on every cell, vectorizes and concanonates

%preallocate
all_heatmapclusters = [];
out_mtx = []; %(4*bins^2, number of cells)
cell_count = 0;

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
            
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_confidence = [3 4 5];
                cluster_region = [1 2];
                cluster_idx = ismember(clusters(:,2), cluster_confidence) & ismember(clusters(:,4), cluster_region);
                clusters = clusters(cluster_idx, :);
                %skip sessions without cells
                if isempty(clusters)
                    continue
                end   
            end 

            %for each cell
            for cf = clusters(:,1)'
                cell_count = cell_count+1;
                rv_cell = nan(bins^2*4,1);

                %for each context
                for visit = sort(context_ids)'
                    visit_num = find(visit == sort(context_ids)');

                    %compute firing rate map
                    %skaggsmoothed, smooth2a, inpaintnans
                    %[~,~,~,rm] = rate_mtx(eptrials(eptrials(:,6)==visit,:), cf, bins); 
                    %rm = inpaint_nans(smooth2a(rm,1));
                    [rm] = trlfree_heatmap(eptrials(eptrials(:,6)==visit,:), cf, bins);
                    rm = inpaint_nans(rm);
                    
                    %vectorize
                    lo_rng = 1+bins^2*(visit_num-1);
                    hi_rng = bins^2*(visit_num);
                    
                    %load (concanonate)
                    rv_cell(lo_rng:hi_rng) = rm(:);
                end
                
                %concan
                out_mtx = [out_mtx rv_cell];

            end
        end
    end
end

%zscore all
out_mtx_z = zscore_mtx(out_mtx);


end