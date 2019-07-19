function [out_bw_first, out_bw_second, out_bw_all] = ALL_pop_pixlecorr(bins)
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

out_bw_first = [];
out_bw_second = [];
out_bw_all = [];

%cell_file_ids_bw = cell(1);
%cell_file_ids_ia = cell(1);

%count_bw = 0;
%count_ia = 0;

%names%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);

%iterate through subjects
for subject = 1:length_subjects

    rat = file_names_subjects{:}(subject,1).name;
    task_folders = {'black_and_white' 'object_arrangement'};
    
    for task = 1:2
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));


        %number of folders
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            %day = session
            
            %load session
            current_file = strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file))
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)))
            
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
           
            %sel_contexts = sort(context_ids); sel_contexts = sel_contexts([1 3 2 4]); %set to abab
            
            %INSERT FUNCTION HERE'
            if task == 1
                [vectors_first, vectors_second, vectors_all] = pop_pixlecorr(eptrials, sort(context_ids), clusters(clusters(:,2)>2,1), bins);

                out_bw_first = cat(2, out_bw_first, vectors_first);
                out_bw_second = cat(2, out_bw_second, vectors_second);
                out_bw_all = cat(2, out_bw_all, vectors_all);
                
            elseif task == 2
                continue
            end

        end
    
    end 
end

end

