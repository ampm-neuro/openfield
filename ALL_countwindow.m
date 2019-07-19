function [all_counts_bw, all_counts_ia, context_bw, context_ia,...
    all_counts_bw_zcr, all_counts_ia_zcr, region_idx_bw, region_idx_ia,...
    all_pos_bw, all_hd_bw, all_pos_ia, all_visit_order_bw,...
    all_visit_order_ia, all_hd_ia, all_obj_ia] = ALL_countwindow(window_duration, varargin)
%find spike counts in windows in each of the contexts (including the 3
%ITI visits)
%
%cells are context visits 1,2,3,4,9,10,11,12 for bw and 5,6,7,8,9,10,11,12 
%for ia
%
%plots pca


% can add varargin to query other columns besides the context id (6)
% e.g., col 5 which is order of presentation
if nargin ==2 && varargin{1} == 5
    col = varargin{1};
    contexts_bw = [1 2 3 4];
    contexts_ia = [5 6 7 8];
elseif nargin ==1
    col = 6;
    contexts_bw = [1 2 3 4 10 11 12];
    contexts_ia = [5 6 7 8 10 11 12];
else
    error('check inputs')
end

iti_time = 2;%min

                
all_counts_bw = cell(length(contexts_bw),1);
all_counts_ia = cell(length(contexts_ia),1);

all_vis_order_bw = cell(length(contexts_bw),1);
all_vis_order_ia = cell(length(contexts_bw),1);

all_xpos_bw = cell(length(contexts_bw),1);
all_ypos_bw = cell(length(contexts_bw),1);
all_hd_bw = cell(length(contexts_bw),1);
all_xpos_ia = cell(length(contexts_ia),1);
all_ypos_ia = cell(length(contexts_ia),1);
all_hd_ia = cell(length(contexts_ia),1);
all_obj_ia = cell(length(contexts_ia),1);

region_idx_bw = [];
region_idx_ia = [];

%get all the things in neurodata folder...
file_list_subjects = dir('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};


window_duration

%iterate through subjects
length_subjects = size(file_names_subjects{:},1);
for subject = 1:length_subjects
    
    %iterate through task folders
    rat = file_names_subjects{:}(subject,1).name
    task_folders = {'black_and_white' 'object_arrangement'};
    for task = 1:2
    
        task
        
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat'));
        strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat')
        
        %iterate through sessions
        length_sessions = size(file_list_sessions,1);
        for session = 1:length_sessions

            session
            
            day_file = file_list_sessions(session,1).name;
            
            %
            context_ids = [];
            
            %load session
            current_file = strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file))
            load(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file)), 'clusters', 'context_ids', 'eptrials')
            
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_confidence = [3 4 5];
                cluster_region = [0 1 2];
                cluster_idx = ismember(clusters(:,2), cluster_confidence) & ismember(clusters(:,4), cluster_region);
                clusters = clusters(cluster_idx, :);
                %skip sessions without cells
                if isempty(clusters)
                    continue
                end   
            end  
            
            %INSERT FUNCTION HERE
            wdws_ctxt = floor(12*60 / window_duration);
            wdws_iti = floor(iti_time / window_duration);
            
            if task == 1
                                
                
                region_idx_bw = [region_idx_bw; clusters(:,4)];
                
                
                %iterate through contexts
                for ctxt_counter = 1:length(contexts_bw)
                    current_context_bw = contexts_bw(ctxt_counter);
                    
                    %calculate spike counts
                    [spike_counts_bw, context_visit_order_hold, xpos_hold_bw, ypos_hold_bw, hd_hold_bw] = slide_rate_window(eptrials(eptrials(:,col)==current_context_bw,:), clusters(clusters(:,2)>=3,1), window_duration);
                    cvo_bw = mode(context_visit_order_hold);

                    %if not iti
                    if ismember(current_context_bw, 1:8)
                        %if not 12min of data, duplicate end
                        if size(spike_counts_bw,1) < wdws_ctxt
                            spike_counts_bw = [spike_counts_bw; repmat(spike_counts_bw(end, :), wdws_ctxt - size(spike_counts_bw,1), 1)];
                            xpos_hold_bw = [xpos_hold_bw; repmat(xpos_hold_bw(end, :), wdws_ctxt - size(xpos_hold_bw,1), 1)];
                            ypos_hold_bw = [ypos_hold_bw; repmat(ypos_hold_bw(end, :), wdws_ctxt - size(ypos_hold_bw,1), 1)];
                            hd_hold_bw = [hd_hold_bw; repmat(hd_hold_bw(end, :), wdws_ctxt - size(hd_hold_bw,1), 1)];
                        end
                        %take only first 12 min
                        spike_counts_bw = spike_counts_bw(1:wdws_ctxt, :);
                        xpos_hold_bw = xpos_hold_bw(1:wdws_ctxt, :);
                        ypos_hold_bw = ypos_hold_bw(1:wdws_ctxt, :);
                        hd_hold_bw = hd_hold_bw(1:wdws_ctxt, :);
                        
                    %if iti    
                    else
                        xpos_hold_bw(:,:) = nan; %no iti positions
                        ypos_hold_bw(:,:) = nan; %no iti positions
                        hd_hold_bw(:,:) = nan; %no iti hd
                        if size(spike_counts_bw,1) < floor(wdws_iti)
                            spike_counts_bw = [spike_counts_bw; repmat(spike_counts_bw(end, :), floor(wdws_iti) - size(spike_counts_bw,1), 1)];
                            xpos_hold_bw = [xpos_hold_bw; repmat(xpos_hold_bw(end, :), floor(wdws_iti) - size(xpos_hold_bw,1), 1)];
                            ypos_hold_bw = [ypos_hold_bw; repmat(ypos_hold_bw(end, :), floor(wdws_iti) - size(ypos_hold_bw,1), 1)];
                            hd_hold_bw = [hd_hold_bw; repmat(hd_hold_bw(end, :), floor(wdws_iti) - size(hd_hold_bw,1), 1)];
                        end
                    end
                    

                    %current number of windows
                    crt_win_ct_bw = size(all_counts_bw{ctxt_counter},1);

                    %take\keep the lowest number of windows seen for any
                    %session\cell in that context
                    if ~isempty(all_counts_bw{ctxt_counter})
                        if size(spike_counts_bw,1) > crt_win_ct_bw
                            spike_counts_bw = spike_counts_bw(1:crt_win_ct_bw,:);
                            xpos_hold_bw = xpos_hold_bw(1:crt_win_ct_bw,:);
                            ypos_hold_bw = ypos_hold_bw(1:crt_win_ct_bw,:);
                            hd_hold_bw = hd_hold_bw(1:crt_win_ct_bw,:);
                        elseif size(spike_counts_bw,1) < crt_win_ct_bw
                            all_counts_bw{ctxt_counter} = all_counts_bw{ctxt_counter}(1:size(spike_counts_bw,1),:);
                            all_vis_order_bw{ctxt_counter} = all_vis_order_bw{ctxt_counter}(1:size(spike_counts_bw,1),:);
                            all_xpos_bw{ctxt_counter} = all_xpos_bw{ctxt_counter}(1:size(xpos_hold_bw,1),:);
                            all_ypos_bw{ctxt_counter} = all_ypos_bw{ctxt_counter}(1:size(ypos_hold_bw,1),:);
                            all_hd_bw{ctxt_counter} = all_hd_bw{ctxt_counter}(1:size(hd_hold_bw,1),:);
                        end
                    end
                    
                    
                    %load into context-appropriate cell
                    all_counts_bw{ctxt_counter} = [all_counts_bw{ctxt_counter} spike_counts_bw];
                    %all_vis_order_bw{ctxt_counter} = [all_vis_order_bw{ctxt_counter} repmat(cvo_bw, size(all_vis_order_bw{ctxt_counter},1), size(spike_counts_bw,2))];
                    all_vis_order_bw{ctxt_counter} = [all_vis_order_bw{ctxt_counter} repmat(cvo_bw, size(spike_counts_bw))];
                    all_xpos_bw{ctxt_counter} = [all_xpos_bw{ctxt_counter} repmat(xpos_hold_bw, 1, size(spike_counts_bw,2))];
                    all_ypos_bw{ctxt_counter} = [all_ypos_bw{ctxt_counter} repmat(ypos_hold_bw, 1, size(spike_counts_bw,2))];
                    all_hd_bw{ctxt_counter} = [all_hd_bw{ctxt_counter} repmat(hd_hold_bw, 1, size(spike_counts_bw,2))];
                    
                end

            elseif task == 2
                %continue
                
                %check if all contexts were visited
                if ~isequal(sort(context_ids), [5 6 7 8]')
                    continue
                end
                
                region_idx_ia = [region_idx_ia; clusters(:,4)];
                
                
                %iterate through contexts
                for ctxt_counter = 1:length(contexts_ia)
                    current_context_ia = contexts_ia(ctxt_counter);
                    
                    %calculate spike counts
                    [spike_counts_ia, context_visit_order_hold_ia, xpos_hold_ia,...
                        ypos_hold_ia, hd_hold_ia_allo, hd_hold_ia_obj]...
                        = slide_rate_window_ia(eptrials(eptrials(:,col)==current_context_ia,:),...
                        clusters(clusters(:,2)>=3,1), window_duration);
                    cvo_ia = mode(context_visit_order_hold_ia);
                    
                    %first = size(spike_counts_ia)
                    %first = size(xpos_hold_ia)
                    
                    %if not iti
                    if ismember(current_context_ia, 1:8)
                        %if not 12min of data, duplicate end
                        if size(spike_counts_ia,1) < wdws_ctxt
                            spike_counts_ia = [spike_counts_ia; repmat(spike_counts_ia(end, :), wdws_ctxt - size(spike_counts_ia,1), 1)];
                            xpos_hold_ia = [xpos_hold_ia; repmat(xpos_hold_ia(end, :), wdws_ctxt - size(xpos_hold_ia,1), 1)];
                            ypos_hold_ia = [ypos_hold_ia; repmat(ypos_hold_ia(end, :), wdws_ctxt - size(ypos_hold_ia,1), 1)];
                            hd_hold_ia_allo = [hd_hold_ia_allo; repmat(hd_hold_ia_allo(end, :), wdws_ctxt - size(hd_hold_ia_allo,1), 1)];
                            hd_hold_ia_obj = [hd_hold_ia_obj; repmat(hd_hold_ia_obj(end, :), wdws_ctxt - size(hd_hold_ia_obj,1), 1)];
                                 %RESUME HERE
                        end
                        %take only first 12 min
                        spike_counts_ia = spike_counts_ia(1:wdws_ctxt, :);
                        xpos_hold_ia = xpos_hold_ia(1:wdws_ctxt, :);
                        ypos_hold_ia = ypos_hold_ia(1:wdws_ctxt, :);
                        hd_hold_ia_allo = hd_hold_ia_allo(1:wdws_ctxt, :);
                        hd_hold_ia_obj = hd_hold_ia_obj(1:wdws_ctxt, :);
                    
                    %if iti    
                    else
                        xpos_hold_ia(:,:) = nan; %no iti positions
                        ypos_hold_ia(:,:) = nan; %no iti positions
                        if size(spike_counts_ia,1) < floor(wdws_iti)
                            spike_counts_ia = [spike_counts_ia; repmat(spike_counts_ia(end, :), floor(wdws_iti) - size(spike_counts_ia,1), 1)];
                            xpos_hold_ia = [xpos_hold_ia; repmat(xpos_hold_ia(end, :), floor(wdws_iti) - size(xpos_hold_ia,1), 1)];
                            ypos_hold_ia = [ypos_hold_ia; repmat(ypos_hold_ia(end, :), floor(wdws_iti) - size(ypos_hold_ia,1), 1)];
                            hd_hold_ia_allo = [hd_hold_ia_allo; repmat(hd_hold_ia_allo(end, :), floor(wdws_iti) - size(hd_hold_ia_allo,1), 1)];
                            hd_hold_ia_obj = [hd_hold_ia_obj; repmat(hd_hold_ia_obj(end, :), floor(wdws_iti) - size(hd_hold_ia_obj,1), 1)];
                        end
                    end
                    
                    %second = size(spike_counts_ia)
                    %second = size(xpos_hold_ia)
                    
                    %current number of windows
                    crt_win_ct_ia = size(all_counts_ia{ctxt_counter},1);

                    %take\keep the lowest number of windows seen for any
                    %session\cell in that context
                    if ~isempty(all_counts_ia{ctxt_counter})
                        if size(spike_counts_ia,1) > crt_win_ct_ia
                            spike_counts_ia = spike_counts_ia(1:crt_win_ct_ia,:);
                            xpos_hold_ia = xpos_hold_ia(1:crt_win_ct_ia,:);
                            ypos_hold_ia = ypos_hold_ia(1:crt_win_ct_ia,:);
                            hd_hold_ia_allo = hd_hold_ia_allo(1:crt_win_ct_ia,:);
                            hd_hold_ia_obj = hd_hold_ia_obj(1:crt_win_ct_ia,:);
                        elseif size(spike_counts_ia,1) < crt_win_ct_ia
                            all_counts_ia{ctxt_counter} = all_counts_ia{ctxt_counter}(1:size(spike_counts_ia,1),:);
                            all_vis_order_ia{ctxt_counter} = all_vis_order_ia{ctxt_counter}(1:size(spike_counts_ia,1),:);
                            all_xpos_ia{ctxt_counter} = all_xpos_ia{ctxt_counter}(1:size(xpos_hold_ia,1),:);
                            all_ypos_ia{ctxt_counter} = all_ypos_ia{ctxt_counter}(1:size(ypos_hold_ia,1),:);
                            all_hd_ia{ctxt_counter} = all_hd_ia{ctxt_counter}(1:size(hd_hold_ia_allo,1),:);
                            all_obj_ia{ctxt_counter} = all_obj_ia{ctxt_counter}(1:size(hd_hold_ia_obj,1),:);
                        end
                    end

                    %third = size(spike_counts_ia)
                    %third = size(xpos_hold_ia)
                
                    
                    %load into context-appropriate cell
                    all_counts_ia{ctxt_counter} = [all_counts_ia{ctxt_counter} spike_counts_ia];
                    all_vis_order_ia{ctxt_counter} = [all_vis_order_ia{ctxt_counter} repmat(cvo_ia, size(spike_counts_ia))];
                    %all_vis_order_ia{ctxt_counter} = [all_vis_order_ia{ctxt_counter} repmat(cvo_ia, size(all_vis_order_ia{ctxt_counter},1), size(spike_counts_ia,2))];
                    all_xpos_ia{ctxt_counter} = [all_xpos_ia{ctxt_counter} repmat(xpos_hold_ia, 1, size(spike_counts_ia,2))];
                    all_ypos_ia{ctxt_counter} = [all_ypos_ia{ctxt_counter} repmat(ypos_hold_ia, 1, size(spike_counts_ia,2))];
                    all_hd_ia{ctxt_counter} = [all_hd_ia{ctxt_counter} repmat(hd_hold_ia_allo, 1, size(spike_counts_ia,2))];
                    all_obj_ia{ctxt_counter} = [all_obj_ia{ctxt_counter} repmat(hd_hold_ia_obj, 1, size(spike_counts_ia,2))];
                end
            end
        end
    end 
end

%make context vectors
context_bw = [];
context_ia = [];
for ci = 1:length(contexts_bw)
    context_bw = [context_bw; repmat(contexts_bw(ci),size(all_counts_bw{ci},1),1)];
    context_ia = [context_ia; repmat(contexts_ia(ci),size(all_counts_ia{ci},1),1)];
end


%reshape to matrices
hold_bw = [];
hold_ia = [];
hold_pos_bw_x = [];
hold_pos_bw_y = [];
hold_hd_bw = [];
hold_pos_ia_x = [];
hold_pos_ia_y = [];

hold_allohd_ia = [];
hold_objvars_ia = [];

hold_cvo_bw = [];
hold_cvo_ia = [];

for celli = 1:length(all_counts_bw)
    hold_bw = [hold_bw; all_counts_bw{celli}];
    hold_ia = [hold_ia; all_counts_ia{celli}];
    hold_pos_bw_x = [hold_pos_bw_x; all_xpos_bw{celli}];
    hold_pos_bw_y = [hold_pos_bw_y; all_ypos_bw{celli}];
    hold_hd_bw = [hold_hd_bw; all_hd_bw{celli}];
    hold_pos_ia_x = [hold_pos_ia_x; all_xpos_ia{celli}];
    hold_pos_ia_y = [hold_pos_ia_y; all_ypos_ia{celli}];
    
    hold_allohd_ia = [hold_allohd_ia; all_hd_ia{celli}];
    hold_objvars_ia = [hold_objvars_ia; all_obj_ia{celli}];
    
    hold_cvo_bw = [hold_cvo_bw; all_vis_order_bw{celli}];
    hold_cvo_ia = [hold_cvo_ia; all_vis_order_ia{celli}];
end

%output
all_counts_bw = hold_bw;
all_pos_bw = cat(3, hold_pos_bw_x, hold_pos_bw_y);
all_counts_ia = hold_ia;
all_hd_bw = hold_hd_bw;
all_pos_ia = cat(3, hold_pos_ia_x, hold_pos_ia_y);
all_visit_order_bw = hold_cvo_bw;
all_visit_order_ia = hold_cvo_ia;
all_hd_ia = hold_allohd_ia;
all_obj_ia = hold_objvars_ia;

%reshape all_obj_ia to set 1 3d pages per cell
all_obj_ia = reshape(all_obj_ia, size(all_obj_ia,1), 8, size(all_obj_ia,2)/8);

%figure function
[all_counts_bw_zcr, all_counts_ia_zcr] = bw_ia_fig(all_counts_bw, all_counts_ia, context_bw, context_ia, contexts_bw, contexts_ia);
    function [all_counts_bw_zcr, all_counts_ia_zcr] = bw_ia_fig(all_counts_bw, all_counts_ia, context_bw, context_ia, contexts_bw, contexts_ia)
        
        %color
        colors = [17 17 193;...
                  38 166 230;...
                  202 16 16;...
                  228 146 40;...
                  50 50 50;...
                  115 115 115;...
                  154 154 154].\255;

        %zscore rates
        all_counts_bw_zcr = all_counts_bw - repmat(mean(all_counts_bw), size(all_counts_bw,1), 1);
        all_counts_bw_zcr = all_counts_bw_zcr.\std(all_counts_bw_zcr);
        all_counts_ia_zcr = all_counts_ia - repmat(mean(all_counts_ia), size(all_counts_ia,1), 1);
        all_counts_ia_zcr = all_counts_ia_zcr.\std(all_counts_ia_zcr);

        
        
        
        %principle components
        [all_counts_bw_zcr_pca, ~, ~, ~, explained_bw_pca] = pca(all_counts_bw_zcr');
        [all_counts_ia_zcr_pca, ~, ~, ~, explained_ia_pca] = pca(all_counts_ia_zcr');
        
        %explained_bw_pca = explained_bw_pca
        %explained_ia_pca = explained_ia_pca

        %plot BW
        %{
        figure; hold on
        count = 0;
        for ctxti = contexts_bw
            count = count+1;
            plot3(all_counts_bw_zcr_pca(context_bw==ctxti,1), all_counts_bw_zcr_pca(context_bw==ctxti,2), all_counts_bw_zcr_pca(context_bw==ctxti,3),'.', 'MarkerSize', 8, 'Color', colors(count,:));
            h_bw{count} = plot3(mean(all_counts_bw_zcr_pca(context_bw==ctxti,1)), mean(all_counts_bw_zcr_pca(context_bw==ctxti,2)), mean(all_counts_bw_zcr_pca(context_bw==ctxti,3)),'.', 'MarkerSize', 60, 'Color', colors(count,:)); 
        end
        xlabel('pc1')
        ylabel('pc2')
        zlabel('pc3')
        title('zscored bw')
        legend show
        legend([h_bw{1}, h_bw{2}, h_bw{3}, h_bw{4}, h_bw{5}, h_bw{6}, h_bw{7}], 'Black1', 'Black2', 'White1', 'White2', 'ITI1', 'ITI2', 'ITI3')
        axis square

        %plot IA
        figure; hold on
        count = 0;
        for ctxti = contexts_ia
            count = count+1;
            plot3(all_counts_ia_zcr_pca(context_ia==ctxti,1), all_counts_ia_zcr_pca(context_ia==ctxti,2), all_counts_ia_zcr_pca(context_ia==ctxti,3),'.', 'MarkerSize', 8, 'Color', colors(count,:)); 
            h_ia{count} = plot3(mean(all_counts_ia_zcr_pca(context_ia==ctxti,1)), mean(all_counts_ia_zcr_pca(context_ia==ctxti,2)), mean(all_counts_ia_zcr_pca(context_ia==ctxti,3)),'.', 'MarkerSize', 60, 'Color', colors(count,:)); 
        end
        xlabel('pc1')
        ylabel('pc2')
        zlabel('pc3')
        title('zscored ia')
        legend show
        legend([h_ia{1}, h_ia{2}, h_ia{3}, h_ia{4}, h_ia{5}, h_ia{6}, h_ia{7}], 'Arng1', 'Arng2', 'Obj1', 'Obj2', 'ITI1', 'ITI2', 'ITI3')
        axis square
        %}

    end






end