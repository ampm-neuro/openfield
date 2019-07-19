function ALL_pos(context)
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

figure; hold on; axis square

%context name index
context_names = {'Black1' 'Black2' 'White1' 'White2' 'Arngmt1' 'Arngmt2' 'Object1' 'Object2'};

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

        %get all the things in subject folder...
        file_list_tasks = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name)));

        %hard coded erasure of irrelevant directory folders
        file_list_tasks(1:2) = [];

        %exclude non-folders
        file_names_tasks = {file_list_tasks([file_list_tasks(:).isdir])};

        %number of folders
        length_tasks = size(file_names_tasks{:},1);

        %iterate through stages %CONTROL WHICH STAGE
        for task = 1:length_tasks

            %print update
            task = file_names_tasks{:}(task,1).name;

            %get all the *.mat in subject folder...
            %strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(task), '/*.mat');
            
            file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(task), '/*.mat'));

            %number of folders
            length_sessions = size(file_list_sessions,1);

            %iterate through sessions
            for session = 1:length_sessions

                day_file = file_list_sessions(session,1).name;
                %day = session

                %load session
                load(strcat('neurodata/',num2str(rat), '/' ,num2str(task), '/', num2str(day_file)));

                colors = [0.0000    0.4470    0.7410; ...
                          0.8500    0.3250    0.0980; ...
                          0.9290    0.6940    0.1250; ...
                          0.4940    0.1840    0.5560; ...
                          0.4660    0.6740    0.1880; ...
                          0.3010    0.7450    0.9330; ...
                          0.6350    0.0780    0.1840];
                
                range = [1.1 0.725 0.46 0.25 .11 .035 0];       
                
                
                if ismember(context, unique(eptrials(:,6)))
                
                    %INSERT FUNCTION HERE
                    plot(eptrials(eptrials(:,6)==context,2), eptrials(eptrials(:,6)==context,3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-')
                    
                    %bullseye
                    %{
                    for i = 1:size(colors,1)
                        
                        
                        plot(eptrials(eptrials(:,6)==context & eptrials(:,10)<range(i), 2), eptrials(eptrials(:,6)==context & eptrials(:,10)<range(i), 3), '.', 'markersize', 1, 'color', colors(i,:))
                    
                    
                    
                    end
                    %}
                    %i = .5;
                    %hold on 
                    %plot(eptrials(eptrials(:,6)==context,2), eptrials(eptrials(:,6)==context,3), '-', 'Color', [.8 .8 .8]);
                    %plot(eptrials(eptrials(:,6)==context & eptrials(:,10)<i, 2), eptrials(eptrials(:,6)==context & eptrials(:,10)<i, 3), '.', 'markersize', 1)

                    
                    
                    title(context_names{context})

                end

            end
        end
    end
end