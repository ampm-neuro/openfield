function [bw_out, ia_out, bw_cells, ia_cells, bw_universe, ia_universe, bw_pop, ia_pop, all_bw_clusters, all_ia_clusters, hold_cell, bw_pop_cell_idx] = ALL_cellcounts
%iterates through all sessions counting cells in each of the possible
%context conditions.

hold_cell = {}; hc_count = 0;

bw_universe = [1 2 3 4;1 3 2 4;1 3 4 2;3 1 2 4;3 1 4 2;3 4 1 2];
ia_universe = flipud(perms(5:8));

bw_cells = zeros(size(bw_universe,1),1);
ia_cells = zeros(size(ia_universe,1),1);

bw_out = cell(size(bw_universe));
    bw_out(ismember(bw_universe,[1 2])) = {'B'};
    bw_out(ismember(bw_universe,[3 4])) = {'W'};
ia_out = cell(size(ia_universe));
    ia_out(ia_universe==5) = {'A1'};
    ia_out(ia_universe==6) = {'A2'};
    ia_out(ia_universe==7) = {'O1'};
    ia_out(ia_universe==8) = {'O2'};
    
bw_pop = 0;
ia_pop = 0;
bw_pop_cell_idx = [];

all_bw_clusters = [];
all_ia_clusters = [];
    

%names%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);

%iterate through subjects
for subject = 1:length_subjects
    
    %print update
    rat = file_names_subjects{:}(subject,1).name;
    
    %task folders
    task_folders = {'black_and_white' 'object_arrangement'};
    
    for task = 1:2
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat'));
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            %day = session
            
            %load session
            %strcat('\Users\ampm\Documents\MATLAB\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file));
            load(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file)));
            
            if isempty(clusters)
                continue
            end
            %{
            %cluster index specifications
            cluster_confidence = [3 4 5];
            %waveform_shape = [0 1 2 3];
            %hemisphere = [0 1];
            cluster_idx = ismember(clusters(:,2), cluster_confidence);% & ismember(clusters(:,3),waveform_shape);
            clusters = clusters(cluster_idx, :);
            %}
                %constrain clusters
                cluster_confidence = [3 4 5];
                cluster_region = [0 1 2];
                cluster_idx = ismember(clusters(:,2), cluster_confidence) & ismember(clusters(:,4), cluster_region);
                clusters = clusters(cluster_idx, :);
                
            %skip sessions without cells
            if isempty(clusters)
                continue
            end
            
            %{
            if isequal([3 4 1 2], context_ids')
                hc_count = hc_count+1
                hold_cell{hc_count} = strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file));
            end
%}
            
            %INSERT FUNCTION HERE'
            %{
            if task == 1
                    bw_cells(ismember(bw_universe, context_ids', 'rows')) = bw_cells(ismember(bw_universe, context_ids', 'rows')) + sum(clusters(:,2)>=3);
            elseif task == 2
                    ia_cells(ismember(ia_universe, context_ids', 'rows')) = ia_cells(ismember(ia_universe, context_ids', 'rows')) + sum(clusters(:,2)>=3);
            end
            %}
            
            if task == 1
                
                
                %check if all contexts were visited
                if ~isequal(sort(context_ids), [1 2 3 4]')
                    continue
                end
                
                    bw_cells(ismember(bw_universe, context_ids', 'rows')) = bw_cells(ismember(bw_universe, context_ids', 'rows')) + size(clusters,1);
                    

                    %clusters
                    all_bw_clusters = [all_bw_clusters; clusters];

                    %population
                    if size(clusters,1) >= 8 %sum(unique(eptrials(~isnan(eptrials(:,4)),4))>1) >= 8
                        bw_pop = bw_pop + 1;
                        bw_pop_cell_idx = [bw_pop_cell_idx; ones(size(clusters(:,1))).*bw_pop];  
                        strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file))
                    else
                        bw_pop_cell_idx = [bw_pop_cell_idx; zeros(size(clusters(:,1)))];
                    end

            elseif task == 2
                
                %check if all contexts were visited
                if ~isequal(sort(context_ids), [5 6 7 8]')
                    continue
                end
                
                    ia_cells(ismember(ia_universe, context_ids', 'rows')) = ia_cells(ismember(ia_universe, context_ids', 'rows')) + size(clusters,1);
            
                    all_ia_clusters = [all_ia_clusters; clusters];
                    
                    %population
                    if sum(unique(eptrials(~isnan(eptrials(:,4)),4))>1) >= 8
                        ia_pop = ia_pop + 1;
                    end
            end
            
            
            
            %clear
            clusters = [];
            context_ids = [];
            
        end
    
    end 
end

%combine character arrays for output
bw_cell_hold = cell(size(bw_out,1), 1);
for ibch = 1:length(bw_cells) 
    bw_cell_hold(ibch) = {num2str(bw_cells(ibch))};
end
bw_out = [bw_out bw_cell_hold];

ia_cell_hold = cell(size(ia_out,1), 1);
for iich = 1:length(ia_cells) 
    ia_cell_hold(iich) = {num2str(ia_cells(iich))};
end
ia_out = [ia_out ia_cell_hold];


end