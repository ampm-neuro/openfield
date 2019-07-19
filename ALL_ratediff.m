function [ratediff_bw, ratediff_ia, cell_file_ids_bw, cell_file_ids_ia] = ALL_ratediff
%find firing rate differences for within and between contexts for BW and
% for within arrangment and within object contexts for IA

ratediff_bw = [];
ratediff_ia = [];

cell_file_ids_bw = cell(1);
cell_file_ids_ia = cell(1);

count_bw = 0;
count_ia = 0;

%names%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');

%hard coded erasure of irrelevant directory folders
file_list_subjects(1:2) = [];

%exclude non-folders
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%number of folders
length_subjects = size(file_names_subjects{:},1);
%file_names_subjects{:}(1:length_subjects,1).name;

%iterate through subjects
for subject = 1:length_subjects
    
    %print update
    rat = file_names_subjects{:}(subject,1).name;
    
    %task folders
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
            
            if isempty(clusters)
                continue
            end            
            
            %INSERT FUNCTION HERE
            if task == 1

                [within, between] = ctxt_rate_diff(eptrials, clusters(:,1), context_ids);    

                ratediff_bw = [ratediff_bw; [within between]];
                
                
                for ci = 1:size(within,1)
                    count_bw = count_bw+1;
                    cell_file_ids_bw{count_bw} = [current_file '/' num2str(clusters(ci,1))];
                end
                

            elseif task == 2
                
                
                
                [within_a, within_o] = ctxt_rate_diff(eptrials, clusters(:,1), context_ids);    

                ratediff_ia = [ratediff_ia; [within_a within_o]];
                
                
                for ci = 1:size(within_a,1)
                    count_ia = count_ia+1;
                    cell_file_ids_ia{count_ia} = [current_file '/' num2str(clusters(ci,1))];
                end
            end

        end
    
    end 
end

cell_file_ids_bw = cell_file_ids_bw';
cell_file_ids_ia = cell_file_ids_ia';

end