function [clusters_bw, clusters_ia, context_ids_bw, context_ids_ia] = ALL_clusters
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

clusters_bw = [];
clusters_ia = [];
context_ids_bw = [];
context_ids_ia = [];

%names%get all the things in neurodata folder...
file_list_subjects = dir('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\');

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
        file_list_sessions = dir(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat'));


        %number of folders
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            %day = session
            
            %load session
            strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file));
            load(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file)))
            
            if isempty(clusters)
                continue
            end
                       
            %INSERT FUNCTION HERE'
            if task == 1
                    clusters_bw = [clusters_bw; clusters];
                    context_ids_bw = [context_ids_bw; repmat(context_ids', size(clusters,1), 1)];

            elseif task == 2
                    clusters_ia = [clusters_ia; clusters];
                    context_ids_ia = [context_ids_ia; repmat(context_ids', size(clusters,1), 1)];
                
            end

        end
    
    end 
end


end

