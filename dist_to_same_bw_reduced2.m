function dists_out = dist_to_same_bw_reduced2(window_duration, aci)
%calculate population distance to same context type (e.g., b to b) as a function of visit
%number


%just keeps count of cells
all_cell_idx = [];


%expected contexts
contexts_bw = [1 2 3 4];% 9 10 11 12 13];

%minimum confidence
min_conf = 3;

%regions
regions = [0 1 2];
%regions = 2;

%preallocate
all_counts_1 = cell(2,1); %first visit, second visit
all_counts_2 = cell(2,1);
all_counts_3 = cell(2,1);

%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%iterate through subjects
length_subjects = size(file_names_subjects{:},1);
for subject = 1:length_subjects
    
    %iterate through task folders
    rat = file_names_subjects{:}(subject,1).name
    task_folders = {'black_and_white' 'object_arrangement'};
    for task = 1 %bw only
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));

        %iterate through sessions
        length_sessions = size(file_list_sessions,1);
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            
            %load session
            current_file = strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file));
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)), 'clusters', 'context_ids', 'eptrials')
            
            low_idx = length(all_cell_idx)+1
            
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_idx = clusters(:,2)>=min_conf & ismember(clusters(:,4), regions);
                clusters = clusters(cluster_idx, :);
                %skip sessions without cells
                if isempty(clusters)
                    continue
                end   
                all_cell_idx = [all_cell_idx; ones(size(clusters(:,1)))];
            end  
            
            %check if all task contexts were visited
            visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
            if ~isequal(visited_contexts(1:4), contexts_bw')
                all_cell_idx = [all_cell_idx; zeros(size(clusters(:,1)))];
                continue
            end
            
            hi_idx = length(all_cell_idx)
            
            
            %only  use  <= 30 clearest cells
            clear_signal_noise_idx = aci(low_idx : hi_idx);
            clear_signal_noise_idx = logical(clear_signal_noise_idx);

            clusters = clusters(clear_signal_noise_idx, :);
            if isempty(clusters)
                    continue
            end 
            

            %get spike counts
            %contexts output as 4 cells in context_ids order
            context_ids = context_ids';
            [spike_count] = get_allcounts(eptrials, clusters, context_ids, window_duration, min_conf);

            %load based on position distance between alike contexts 
            %(WITHIN CONTEXTS)
            %
            switch 1
                case ismember(context_ids, [1 2 3 4; 3 4 1 2], 'rows') %AABB
                    
  
                    [all_counts_1{1}, spike_count{1}] = uniformtw(all_counts_1{1}, spike_count{1});
                    [all_counts_1{2}, spike_count{2}] = uniformtw(all_counts_1{2}, spike_count{2});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([1 2]));
                                        
                    [all_counts_1{1}, spike_count{3}] = uniformtw(all_counts_1{1}, spike_count{3});
                    [all_counts_1{2}, spike_count{4}] = uniformtw(all_counts_1{2}, spike_count{4});    
                    all_counts_1 = conc_cells(all_counts_1, spike_count([3 4]));
                
                case ismember(context_ids, [1 3 2 4; 3 1 4 2], 'rows') %ABAB
                    
                    [all_counts_2{1}, spike_count{1}] = uniformtw(all_counts_2{1}, spike_count{1});
                    [all_counts_2{2}, spike_count{3}] = uniformtw(all_counts_2{2}, spike_count{3}); 
                    all_counts_2 = conc_cells(all_counts_2, spike_count([1 3]));
                    
                    [all_counts_2{1}, spike_count{2}] = uniformtw(all_counts_2{1}, spike_count{2});
                    [all_counts_2{2}, spike_count{4}] = uniformtw(all_counts_2{2}, spike_count{4});                     
                    all_counts_2 = conc_cells(all_counts_2, spike_count([2 4]));
                
                case ismember(context_ids, [1 3 4 2; 3 1 2 4], 'rows') %ABBA
                    
                    [all_counts_1{1}, spike_count{2}] = uniformtw(all_counts_1{1}, spike_count{2});
                    [all_counts_1{2}, spike_count{3}] = uniformtw(all_counts_1{2}, spike_count{3});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([2 3]));
                    
                    [all_counts_3{1}, spike_count{1}] = uniformtw(all_counts_3{1}, spike_count{1});
                    [all_counts_3{2}, spike_count{4}] = uniformtw(all_counts_3{2}, spike_count{4});
                    all_counts_3 = conc_cells(all_counts_3, spike_count([1 4]));
            end
            %}
            
            %load based on position distance between distinct contexts
            %(BETWEEN CONTEXTS)
            %{
            switch 1
                case ismember(context_ids, [1 2 3 4; 3 4 1 2], 'rows') %AABB
                    
  
                    [all_counts_1{1}, spike_count{2}] = uniformtw(all_counts_1{1}, spike_count{2});
                    [all_counts_1{2}, spike_count{3}] = uniformtw(all_counts_1{2}, spike_count{3});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([2 3]));
                                        
                    [all_counts_2{1}, spike_count{1}] = uniformtw(all_counts_2{1}, spike_count{1});
                    [all_counts_2{2}, spike_count{3}] = uniformtw(all_counts_2{2}, spike_count{3}); 
                    all_counts_2 = conc_cells(all_counts_2, spike_count([1 3]));
                    
                    [all_counts_2{1}, spike_count{2}] = uniformtw(all_counts_2{1}, spike_count{2});
                    [all_counts_2{2}, spike_count{4}] = uniformtw(all_counts_2{2}, spike_count{4});                     
                    all_counts_2 = conc_cells(all_counts_2, spike_count([2 4]));
                    
                    [all_counts_3{1}, spike_count{1}] = uniformtw(all_counts_3{1}, spike_count{1});
                    [all_counts_3{2}, spike_count{4}] = uniformtw(all_counts_3{2}, spike_count{4});
                    all_counts_3 = conc_cells(all_counts_3, spike_count([1 4]));
                    
                    
                    
                case ismember(context_ids, [1 3 2 4; 3 1 4 2], 'rows') %ABAB
                    
                    [all_counts_1{1}, spike_count{1}] = uniformtw(all_counts_1{1}, spike_count{1});
                    [all_counts_1{2}, spike_count{2}] = uniformtw(all_counts_1{2}, spike_count{2});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([1 2]));
                                        
                    [all_counts_1{1}, spike_count{3}] = uniformtw(all_counts_1{1}, spike_count{3});
                    [all_counts_1{2}, spike_count{4}] = uniformtw(all_counts_1{2}, spike_count{4});    
                    all_counts_1 = conc_cells(all_counts_1, spike_count([3 4]));
                    
                    [all_counts_1{1}, spike_count{2}] = uniformtw(all_counts_1{1}, spike_count{2});
                    [all_counts_1{2}, spike_count{3}] = uniformtw(all_counts_1{2}, spike_count{3});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([2 3]));
                    
                    [all_counts_3{1}, spike_count{1}] = uniformtw(all_counts_3{1}, spike_count{1});
                    [all_counts_3{2}, spike_count{4}] = uniformtw(all_counts_3{2}, spike_count{4});
                    all_counts_3 = conc_cells(all_counts_3, spike_count([1 4]));
                
                case ismember(context_ids, [1 3 4 2; 3 1 2 4], 'rows') %ABBA
                    
                    [all_counts_1{1}, spike_count{1}] = uniformtw(all_counts_1{1}, spike_count{1});
                    [all_counts_1{2}, spike_count{2}] = uniformtw(all_counts_1{2}, spike_count{2});
                    all_counts_1 = conc_cells(all_counts_1, spike_count([1 2]));
                                        
                    [all_counts_1{1}, spike_count{3}] = uniformtw(all_counts_1{1}, spike_count{3});
                    [all_counts_1{2}, spike_count{4}] = uniformtw(all_counts_1{2}, spike_count{4});    
                    all_counts_1 = conc_cells(all_counts_1, spike_count([3 4]));
                    
                    [all_counts_2{1}, spike_count{1}] = uniformtw(all_counts_2{1}, spike_count{1});
                    [all_counts_2{2}, spike_count{3}] = uniformtw(all_counts_2{2}, spike_count{3}); 
                    all_counts_2 = conc_cells(all_counts_2, spike_count([1 3]));
                    
                    [all_counts_2{1}, spike_count{2}] = uniformtw(all_counts_2{1}, spike_count{2});
                    [all_counts_2{2}, spike_count{4}] = uniformtw(all_counts_2{2}, spike_count{4});                     
                    all_counts_2 = conc_cells(all_counts_2, spike_count([2 4]));
                    
            end
            %}
            
        end
    end 
end

%COMPUTE DISTANCES
%
dists_out{1} = popdist(all_counts_1);
dists_out{2} = popdist(all_counts_2);
dists_out{3} = popdist(all_counts_3);


%PLOT
%
errorbar_barplot( dists_out )
xlabel('Distance between positions')
ylabel('Population distance')


%INTERNAL FUNCTIONS
%
%allcounts function
    function [all_counts] = get_allcounts(eptrials, clusters, task_context_ids, window_duration, min_conf)

        %preallocate
        all_counts = cell(length(task_context_ids), 1);
        
        %calculate spike counts
        [spike_counts, context_idx] = slide_rate_window(eptrials, clusters(clusters(:,2)>=min_conf, 1), window_duration);
        
        %zscore spike counts
        spike_counts = zscore_mtx(spike_counts);

        %iterate through contexts
        for ctxt_counter = 1:length(task_context_ids)
            current_context = task_context_ids(ctxt_counter);

            %load into context-appropriate cell
            all_counts{ctxt_counter} = [all_counts{ctxt_counter} spike_counts(context_idx==current_context, :)];
        end
    end

%elementwise cocantenation of cells
    function C = conc_cells(A, B)
        C = cell(size(A));
        for cl = 1:numel(A)
            C{cl} = [A{cl} B{cl}];
        end
    end

%pop distance function
%
    function dists_out = popdist(clouds)
        %calculate distances from all samples to mean of other cloud
        
        %preallocate output
        dists_out = [];
        
        %iterate through clouds
        for cloud = 1:length(clouds)
            
            %oposite mean
            if cloud == 1
                opposite_mean = mean(clouds{2});
            else
                opposite_mean = mean(clouds{1});
            end

            %calculate dist from every samp to opposite mean
            dists_to_mean = pdist([opposite_mean; clouds{cloud}]); 
            dists_to_mean = dists_to_mean(1:size(clouds{cloud},1));%ignore extraneous
            dists_to_mean = dists_to_mean./sqrt(size(clouds{cloud},2)); %correct for dimensionality
            dists_to_mean = dists_to_mean - 1; %reset y axis to pos/neg around 1

            %load output
            dists_out = [dists_out; dists_to_mean'];
        end
    end

%uniform window lengths
%
    function [exist_out, new_out] = uniformtw(exist_mtx, new_mtx)
        %sets both input matrices to have the same number of rows as the
        %matrix with the least number of rows
        
        %preemptively set
        exist_out = exist_mtx;
        new_out = new_mtx;
        
        %current number of windows
        crt_win_ct_bw = size(exist_mtx,1);

        %take/keep the lowest number of windows seen for any
        %session/cell in that context
        if ~isempty(exist_mtx)
            if size(new_mtx,1) > crt_win_ct_bw
                new_out = new_mtx(1:crt_win_ct_bw,:);
            elseif size(new_mtx,1) < crt_win_ct_bw
                exist_out = exist_mtx(1:size(new_mtx,1),:);
            end
        end
        
    end

end