function [rm_bw, rm_ia, rm_bw_z, rm_ia_z, context_bw, context_ia, first_second_idx] = ALL_countwindow_prepost(window_duration, context_or_visit)
%find spike counts in windows in each of the contexts (including the 3
%ITI visits)
%
%cells are context visits 1,2,3,4,9,10,11,12 for bw and 5,6,7,8,9,10,11,12 
%for ia
%
%plots pca
%
%context_or_visit
%calculate spike count matrices grouped by either context or visit number
%context is 1, visit is 2


% can add varargin to query other columns besides the context id (6)
% e.g., col 5 which is order of presentation
contexts_bw = [1 2 3 4];% 9 10 11 12 13];
contexts_ia = [5 6 7 8];% 9 10 11 12 13];



%time included from itis around context visits 
pre_post_time = 60;%seconds

%preallocate
all_counts_bw = cell(length(contexts_bw), 2);
all_counts_ia = cell(length(contexts_ia), 2);

%times to keep
before_after_time = 30; %seconds
visit_time = 6*60;

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
    for task = 1:2
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));

        %iterate through sessions
        length_sessions = size(file_list_sessions,1);
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            
            %
            context_ids = [];
            
            %load session
            current_file = strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file));
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)), 'clusters', 'context_ids', 'eptrials')
            
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_confidence = [3 4 5];
                cluster_idx = ismember(clusters(:,2), cluster_confidence);
                clusters = clusters(cluster_idx, :);
                %skip sessions without cells
                if isempty(clusters)
                    continue
                end   
            end  
            
            %skip sessions with too short beginings or ends
            if (min(eptrials(eptrials(:,5)==1,1))-pre_post_time) < 0 | (max(eptrials(eptrials(:,5)==4,1))+pre_post_time) > max(eptrials(:,1)) %#ok<NODEF,OR2>
                continue
            end

            %get spike counts
            if task == 1
                     
                %check if all task contexts were visited
                visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
                if ~isequal(visited_contexts(1:4), contexts_bw')
                    continue
                end
                
                all_counts_bw = conc_cells(all_counts_bw, get_allcounts(eptrials, clusters, contexts_bw, window_duration, pre_post_time, context_or_visit));


            elseif task == 2
                
                %check if all task contexts were visited
                visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
                if ~isequal(visited_contexts(1:4), contexts_ia')
                    continue
                end

                all_counts_ia = conc_cells(all_counts_ia, get_allcounts(eptrials, clusters, contexts_ia, window_duration, pre_post_time, context_or_visit));
            end
        end
    end 
end


%preallocate ids
context_bw = [];
context_ia = [];
first_second_idx = [];
rm_bw = [];
rm_ia = [];

%firsthalf_2ndhalf
for fhsh = 1:2
    
    %make context indices
    context_bw = [];
    context_ia = [];

    for ci = 1:length(contexts_bw)
        context_bw = [context_bw; repmat([contexts_bw(ci) ci], size(all_counts_bw{ci, fhsh}, 1), 1)];
        context_ia = [context_ia; repmat([contexts_ia(ci) ci], size(all_counts_ia{ci, fhsh}, 1), 1)];
    end

    %reshape cells to matrices
    for celli = 1:size(all_counts_bw,1)
        rm_bw = [rm_bw; all_counts_bw{celli, fhsh}];
        rm_ia = [rm_ia; all_counts_ia{celli, fhsh}];
    end
    
    %firsthalf_2ndhalf idx
    first_second_idx = [first_second_idx; repmat(fhsh,size(context_bw))];

end


%zscore rates
rm_bw_z = rm_bw - repmat(mean(rm_bw), size(rm_bw,1), 1);
rm_bw_z = rm_bw_z./std(rm_bw_z);
rm_ia_z = rm_ia - repmat(mean(rm_ia), size(rm_ia,1), 1);
rm_ia_z = rm_ia_z./std(rm_ia_z);

%PLOT FIGURE
context_bw = context_bw(:,context_or_visit); 
context_ia = context_ia(:,context_or_visit);
bw_ia_fig(rm_bw_z, rm_ia_z, context_bw, context_ia, first_second_idx);



%INTERNAL FUNCTIONS
%
%allcounts function
    function all_counts = get_allcounts(eptrials, clusters, task_context_ids, window_duration, pre_post_time, ctx_or_vis)

        %preallocate
        all_counts = cell(length(task_context_ids), 2);

        %befores then afters
        for firsthalf_2ndhalf = 1:2

            %iterate through contexts
            for ctxt_counter = 1:length(task_context_ids)
                current_context = task_context_ids(ctxt_counter);

                %key time
                %
                if ctx_or_vis == 1
                    %Context
                    if firsthalf_2ndhalf ==1
                        keytime = min(eptrials(eptrials(:,6)==current_context, 1));
                        keytime_start = keytime - pre_post_time;
                        keytime_end = keytime+(6*60);
                    elseif firsthalf_2ndhalf ==2
                        keytime = max(eptrials(eptrials(:,6)==current_context, 1));
                        keytime_start = keytime -(6*60);  
                        keytime_end = keytime + pre_post_time;
                    end
                elseif ctx_or_vis == 2
                    %Visit number
                    if firsthalf_2ndhalf ==1
                        keytime = min(eptrials(eptrials(:,5)==ctxt_counter, 1));
                        keytime_start = keytime - pre_post_time;
                        keytime_end = keytime+(6*60);
                    elseif firsthalf_2ndhalf ==2
                        keytime = max(eptrials(eptrials(:,5)==ctxt_counter, 1));
                        keytime_start = keytime -(6*60);  
                        keytime_end = keytime + pre_post_time;
                    end
                end
                keytime_idx = eptrials(:,1) > keytime_start  &  eptrials(:,1) < keytime_end;
                
                %calculate spike counts
                [spike_counts] = slide_rate_window(eptrials(keytime_idx,:), clusters(clusters(:,2)>=3,1), window_duration);

                %load into context-appropriate cell
                all_counts{ctxt_counter, firsthalf_2ndhalf} = [all_counts{ctxt_counter, firsthalf_2ndhalf} spike_counts];
            end
        end
    end

%elementwise cocantenation of cells
    function C = conc_cells(A, B)
        C = cell(size(A));
        for cl = 1:numel(A)
            C{cl} = [A{cl} B{cl}];
        end
    end


%figure function
%
    function bw_ia_fig(rm_bw, rm_ia, ctx_idx_bw, ctx_idx_ia, first_second_idx)
        
        %color
        colors = [17 17 193;...
                  38 166 230;...
                  202 16 16;...
                  228 146 40;...
                  50 50 50;...
                  115 115 115;...
                  154 154 154]./255;
        
        %principle components
        [rm_bw_pca, ~, ~, ~, explained_bw_pca] = pca(rm_bw');
        [rm_ia_pca, ~, ~, ~, explained_ia_pca] = pca(rm_ia');
        
        explained_bw_pca = sum(explained_bw_pca(1:3))
        explained_ia_pca = sum(explained_ia_pca(1:3))

        %plot BW
        figure; hold on
        count = 0;
        for ctxti = unique(ctx_idx_bw)'
            count = count+1;
            plot3(rm_bw_pca(ctx_idx_bw==ctxti,1), rm_bw_pca(ctx_idx_bw==ctxti,2), rm_bw_pca(ctx_idx_bw==ctxti,3), '.', 'MarkerSize', 8, 'Color', colors(count,:));
            h_bw{count} = plot3(mean(rm_bw_pca(ctx_idx_bw==ctxti,1)), mean(rm_bw_pca(ctx_idx_bw==ctxti,2)), mean(rm_bw_pca(ctx_idx_bw==ctxti,3)), '.', 'MarkerSize', 60, 'Color', colors(count,:)); 
        end
        xlabel('pc1')
        ylabel('pc2')
        zlabel('pc3')
        title('zscored bw')
        legend show
        %legend([h_bw{1}, h_bw{2}, h_bw{3}, h_bw{4}, h_bw{5}, h_bw{6}, h_bw{7}], 'Black1', 'Black2', 'White1', 'White2', 'ITI1', 'ITI2', 'ITI3')
        %legend([h_bw{1}, h_bw{2}, h_bw{3}, h_bw{4}], 'Black1', 'Black2', 'White1', 'White2')
        legend([h_bw{1}, h_bw{2}, h_bw{3}, h_bw{4}], 'Visit1', 'Visit2', 'Visit3', 'Visit4')

        axis square

        %plot IA
        figure; hold on
        count = 0;
        for ctxti = unique(ctx_idx_ia)'
            count = count+1;
            plot3(rm_ia_pca(ctx_idx_ia==ctxti,1), rm_ia_pca(ctx_idx_ia==ctxti,2), rm_ia_pca(ctx_idx_ia==ctxti,3), '.', 'MarkerSize', 8, 'Color', colors(count,:)); 
            h_ia{count} = plot3(mean(rm_ia_pca(ctx_idx_ia==ctxti,1)), mean(rm_ia_pca(ctx_idx_ia==ctxti,2)), mean(rm_ia_pca(ctx_idx_ia==ctxti,3)), '.', 'MarkerSize', 60, 'Color', colors(count,:)); 
        end
        xlabel('pc1')
        ylabel('pc2')
        zlabel('pc3')
        title('zscored ia')
        legend show
        %legend([h_ia{1}, h_ia{2}, h_ia{3}, h_ia{4}, h_ia{5}, h_ia{6}, h_ia{7}], 'Arng1', 'Arng2', 'Obj1', 'Obj2', 'ITI1', 'ITI2', 'ITI3')
        legend([h_ia{1}, h_ia{2}, h_ia{3}, h_ia{4}], 'Visit1', 'Visit2', 'Visit3', 'Visit4')
        axis square

    end





end