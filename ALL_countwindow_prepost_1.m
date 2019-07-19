function [rm, rm_z, ctxid, visit_mtx, first_second_idx] = ALL_countwindow_prepost_1(window_duration, pre_time, iti_time, post_time, aci)
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

all_cell_idx = [];


% can add varargin to query other columns besides the context id (6)
% e.g., col 5 which is order of presentation
contexts_bw = [1 2 3 4];% 9 10 11 12 13];
contexts_ia = [5 6 7 8];% 9 10 11 12 13];

%first tw start buffer
strt_buf = 10;

%minimum confidence
min_conf = 0;

%regions
regions = [0 1 2];
%regions = 1;

%visit time
visit_time = 6*60;%seconds
iti_time = floor(iti_time/2);

%time included from itis around context visits %
%pre_post_time = 60;%seconds

%preallocate
all_counts = cell(length(contexts_bw), 2);
all_ctx_id = cell(length(contexts_bw), 2);
all_first_tw = [];

%get all the things in neurodata folder...
file_list_subjects = dir('neurodata\');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%iterate through subjects
length_subjects = size(file_names_subjects{:},1);
for subject = 1:length_subjects
    
    %iterate through task folders
    rat = file_names_subjects{:}(subject,1).name;
    task_folders = {'black_and_white' 'object_arrangement'};
    for task = 1:2
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat'));

        %iterate through sessions
        length_sessions = size(file_list_sessions,1);
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            
            %load session
            current_file = strcat('\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file));
            load(strcat('\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file)), 'clusters', 'context_ids', 'eptrials')
                        
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_confidence = [3 4 5];
                cluster_region = [0 1 2];
                cluster_idx = ismember(clusters(:,2), cluster_confidence) & ismember(clusters(:,4), cluster_region);
                clusters = clusters(cluster_idx, :);
                if isempty(clusters)
                    continue
                end 
            end  
            
            low_idx = length(all_cell_idx)+1;
            
            %skip sessions with too short beginings or ends
            if (min(eptrials(eptrials(:,5)==1,1))-pre_time) < 0 | (max(eptrials(eptrials(:,5)==4,1))+post_time) > max(eptrials(:,1)) %#ok<NODEF,OR2>
                all_cell_idx = [all_cell_idx; zeros(size(clusters(:,1)))];
                continue
            else
                all_cell_idx = [all_cell_idx; ones(size(clusters(:,1)))];
            end
            
            hi_idx = length(all_cell_idx);

            %only  use  30 clearest cells (24 after other screens)
            %{
            
            clear_signal_noise_idx = aci(low_idx : hi_idx);
            clear_signal_noise_idx = logical(clear_signal_noise_idx);

            clusters = clusters(clear_signal_noise_idx, :);
            if isempty(clusters)
                    continue
            end 
            %}
            
            
            %get spike counts
            if task == 1
                     
                %check if all task contexts were visited
                visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
                if ~isequal(visited_contexts(1:4), contexts_bw')
                    continue
                end

                %first sample
                first_tw = slide_rate_window(eptrials(eptrials(:,1)>strt_buf & eptrials(:,1)<=strt_buf+window_duration+1,:), clusters(clusters(:,2)>=min_conf,1), window_duration);
                first_tw = first_tw(1,:);
                
                %full session
                [spike_count, tw_id] = get_allcounts(eptrials, clusters, contexts_bw, window_duration, pre_time, iti_time, post_time, min_conf);
                
                %load
                all_first_tw = [all_first_tw first_tw];
                all_counts = conc_cells(all_counts, spike_count);
                all_ctx_id = conc_cells(all_ctx_id, tw_id);


            elseif task == 2
                
                %IF JUST BW, USE THIS CONTINUE
                %continue
                
                %check if all task contexts were visited
                visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
                if ~isequal(visited_contexts(1:4), contexts_ia')
                    continue
                end
                
                %first sample
                first_tw = slide_rate_window(eptrials(eptrials(:,1)>strt_buf & eptrials(:,1)<=strt_buf+window_duration+1,:), clusters(clusters(:,2)>=min_conf,1), window_duration);
                first_tw = first_tw(1,:);
                
                %full session
                [spike_count, tw_id] = get_allcounts(eptrials, clusters, contexts_ia, window_duration, pre_time, iti_time, post_time, min_conf);
                
                %load
                all_first_tw = [all_first_tw first_tw];
                all_counts = conc_cells(all_counts, spike_count);
                all_ctx_id = conc_cells(all_ctx_id, tw_id);
            end
            %}
        end
    end 
end



%preallocate ids
visit_idx = [];
first_second_idx = [];
rm = [];
ctxid = [];

%firsthalf_2ndhalf
for fhsh = 1:2

    for ci = 1:4
        visit_idx = [visit_idx; repmat(ci, size(all_counts{ci, fhsh}, 1), 1)];
    end

    %reshape spike count cells to matrices
    rm_hold = [];
    for celli = 1:size(all_counts,1)
        rm_hold = [rm_hold; all_counts{celli, fhsh}];
    end
    rm = [rm; rm_hold];
    
    
    %reshape context id cells to matrices
    ctxid_hold = [];
    for celli = 1:size(all_ctx_id,1)
        ctxid_hold = [ctxid_hold; all_ctx_id{celli, fhsh}];
    end
    ctxid = [ctxid; ctxid_hold];
    
    %firsthalf_2ndhalf idx
    first_second_idx = [first_second_idx; repmat(fhsh, size(rm_hold))];

end

%failed idea (remains throughout code)
all_first_tw = [];
first_tw_z = [];

%zscore rates
rm_hold = [all_first_tw; rm];
rm_z = rm_hold - repmat(mean(rm_hold), size(rm_hold,1), 1);
rm_z = rm_z./std(rm_z);

%first_tw_z = rm_z(1:size(all_first_tw,1),:);
%rm_z = rm_z((size(all_first_tw,1)+1) : end,:);

%sort rm outputs
rm_z_hold = [];
rm_hold = [];
ctxid_hold = [];
first_second_idx_hold = [];
visit_idx_hold = [];
for vi = unique(visit_idx)'
    rm_hold = [rm_hold; rm(visit_idx==vi, :)];
    rm_z_hold = [rm_z_hold; rm_z(visit_idx==vi, :)];
    ctxid_hold = [ctxid_hold; ctxid(visit_idx==vi, :)];
    first_second_idx_hold = [first_second_idx_hold; first_second_idx(visit_idx==vi, :)];
    visit_idx_hold = [visit_idx_hold; visit_idx(visit_idx==vi, :)];
end
rm = rm_hold;
rm_z = rm_z_hold;
ctxid = ctxid_hold;
first_second_idx = first_second_idx_hold;
visit_idx = visit_idx_hold;

%IDENTIFY visit matrix (when enter and exit contexts)
visit_mtx = [(1:length(mode(ctxid,2)))', mode(ctxid,2)];
visit_mtx = redlines(visit_mtx, 0); %1 for plot on, 0 for off


%PLOT PCA FIGURE
%
comb_fig(rm_z, visit_idx, first_second_idx, visit_mtx); a1 = axis;
figure;
delay_pca_plot(rm_z, 1, 1:3, 51); axis(a1);
%}


%PLOT DRIFT FIGURES
%

    %first_tw of first context    
    first_tw_z = mean(rm_z(find(visit_mtx(:,3)==2, 60/window_duration,'first'), :));

    %dist between all timewindows
    pop_dists = pdist([first_tw_z; rm_z])./sqrt(size(rm_z,2));
    pop_dists = pop_dists(1:size(rm_z,1));
    
    figure; plot(pop_dists); title('popdist test')

    lpd = length(pop_dists);
    %visit_mtx = [(1:length(mode(ctxid,2)))', mode(ctxid,2)];

    figure; hold on
    set(gca,'TickLength',[0, 0]); box off
    
    %red lines marking visit start and end
    %visit_mtx = redlines(visit_mtx, 1); %1 for plot on, 0 for off
    pd_constrained = nan(size(pop_dists,1));
    pd_constrained_prm = nan(size(pop_dists,1));
    pd_constrained(visit_mtx(:,3)==2) = pop_dists(visit_mtx(:,3)==2);
    pd_constrained_prm(visit_mtx(:,3)~=2) = pop_dists(visit_mtx(:,3)~=2);
    
    %line 
    %{
    for il = 2:(length(pop_dists)-1)
        plot((il-1):il, pop_dists(il:il+1,1), 'k-', 'linewidth', 1)
    end
    %}
    
    %dots
    %plot(pd_constrained_prm, 'ko', 'markersize', 10/pi)%, 'Color', [.85 .85 .85])
    plot(pd_constrained, 'k.', 'markersize', 10)
    title('Drift over session')

    xlim([-lpd*.1 lpd*1.1])
    total_time = pre_time + 8*visit_time + 6*iti_time + post_time; %double is normal
    xticks(1:round(lpd/12):lpd)
    xticklabels(round((((2:round(lpd/12):lpd)./lpd).*total_time)./60, 1))
    xlabel('approximate time from start (minutes)')

    ylim([min(pd_constrained(~isnan(pd_constrained)))*.95 max(pd_constrained(~isnan(pd_constrained)))*1.05])
    ylabel('population distance from first time window')

    %fitline and correllations
    avg_ctx = [];
    for cv = unique(visit_mtx(~isnan(visit_mtx(:,4)),4))'
        pd_constrained_fl = nan(size(pop_dists,1));
        pd_constrained_fl(visit_mtx(:,4)==cv) = pop_dists(visit_mtx(:,4)==cv);
        [cv_Rp(cv,1), cv_Rp(cv,2)] = fit_line((1:length(pd_constrained_fl))', pd_constrained_fl);
        
        
        %ERROR HERE Dimensions of arrays being concatenated are not consistent.
        
        
        size((1:sum(~isnan(pd_constrained_fl)))')
        size(zscore_mtx(pd_constrained_fl(~isnan(pd_constrained_fl))))
        

        size([(1:sum(~isnan(pd_constrained_fl)))' zscore_mtx(pd_constrained_fl(~isnan(pd_constrained_fl)))])
        
        avg_ctx = [avg_ctx; [(1:sum(~isnan(pd_constrained_fl)))' zscore_mtx(pd_constrained_fl(~isnan(pd_constrained_fl)))]];
        
    end

    figure; hold on; [all_R, all_p] = fit_line((1:length(pd_constrained))',pd_constrained);
    
    %within session average drift and correllation
    figure; plot(avg_ctx(:,1), avg_ctx(:,2), 'k.')
    [avg_R, avg_p] = fit_line(avg_ctx(:,1), avg_ctx(:,2));
    set(gca,'TickLength',[0, 0]); box off
    title('Drift over visit (average)')
%}
    

    
    
%PLOT POP VELOCITY FIGURE At entrance to visit
%
    distmtx = dist([first_tw_z; rm_z]')./sqrt(size(rm_z,2));
    
    distmtx_mask = diag_mask([size(distmtx,1)+2, size(distmtx,2)+2]); 
    distmtx_mask = distmtx_mask(1:end-2, 3:end); 
    dveloc = nan(1, size([first_tw_z; rm_z],1));
    dveloc(2:(end-1)) = distmtx(distmtx_mask)./(2*window_duration); % /2*window duration = velocity
    [~, ent_ext] = redlines(visit_mtx, 0);
    tw_hilo = [30 60*2];
    dveloc_avg = [];
    for ivdv = 1:4
        dveloc_avg = [dveloc_avg; dveloc(ent_ext(ivdv,1)-tw_hilo(1) : ent_ext(ivdv,1)+tw_hilo(2)-1)];
    end    

    figure; hold on;
    set(gca,'TickLength',[0, 0]); box off
    title('Drift velocity entering visit (average)')
    dveloc_mean = nanmean(dveloc_avg);
    plot(smooth(dveloc_mean,12), 'k-', 'linewidth', 2)
    pos_std = smooth(dveloc_mean+(nanstd(dveloc_avg)./sqrt(4)),12);
    neg_std = smooth(dveloc_mean-(nanstd(dveloc_avg)./sqrt(4)),12);
    plot(pos_std, 'Color', [.5 .5 .5], 'linewidth', 1)
    plot(neg_std, 'Color', [.5 .5 .5], 'linewidth', 1)
    xlim([0 length(dveloc_avg)+1])
    xticks([1 10:10:sum(tw_hilo)])
    xlabel('time(s)')
    ylim([min(neg_std)*.96 max(pos_std)*1.04])
    plot([tw_hilo(1) tw_hilo(1)], ylim, 'r-', 'linewidth', 2)
%}
    
%PLOT DIST TO MEAN FIGURE At entrance to visit
%
    tw_hilo = [30 60*2];
    d2m_avg = [];
    for ivdv = 1:4

        visit_mean = mean(rm_z(visit_mtx(:,4)==ivdv, :));
        rlvt_rates = rm_z(ent_ext(ivdv,1)-tw_hilo(1) : ent_ext(ivdv,1)+tw_hilo(2)-1, :);
        dist_vect = pdist([visit_mean; rlvt_rates])./sqrt(size(rm_z,2));
        dist_vect = dist_vect(1:size(rlvt_rates,1));
        d2m_avg = [d2m_avg; dist_vect];
    end    

    figure; hold on;
    set(gca,'TickLength',[0, 0]); box off
    title('Distance to Visit Mean, entering visit (average)')
    d2m_mean = nanmean(d2m_avg);
    plot(smooth(d2m_mean,12), 'k-', 'linewidth', 2)
    pos_std = smooth(d2m_mean+(nanstd(d2m_avg)./sqrt(4)),12);
    neg_std = smooth(d2m_mean-(nanstd(d2m_avg)./sqrt(4)),12);
    plot(pos_std, 'Color', [.5 .5 .5], 'linewidth', 1)
    plot(neg_std, 'Color', [.5 .5 .5], 'linewidth', 1)
    xlim([0 length(d2m_avg)+1])
    xticks([1 10:10:sum(tw_hilo)])
    xlabel('time(s)')
    ylim([min(neg_std)*.96 max(pos_std)*1.04])
    plot([tw_hilo(1) tw_hilo(1)], ylim, 'r-', 'linewidth', 2)
%}    
    
    
%Correllation of distance-to-mean points upon entrace
% UNINTERESTING
%{
box_car_sz = 30; %points

d2m_mean = smooth(d2m_mean,10);
d2m_mean = smooth(d2m_mean,10);
figure; plot(d2m_mean); title smoothedmean

corr_outputs = nan(length(d2m_mean(tw_hilo(1)+1 : end-box_car_sz)), 2);

for ictw = 1:size(corr_outputs,1)
    crng = tw_hilo(1)+ictw:tw_hilo(1)+ictw+box_car_sz;
    [cR, cp] = corr(crng', d2m_mean(crng));
    corr_outputs(ictw,:) = [cR cp];
end

figure
plot(corr_outputs(:,1), 'k-')
title('Rvalues')
set(gca,'TickLength',[0, 0]); box off

figure 
plot(corr_outputs(:,2), 'k-')
title('pvalues')
set(gca,'TickLength',[0, 0]); box off
%}


num_qualified_cells = length(all_cell_idx)

%INTERNAL FUNCTIONS
%
%allcounts function
    function [all_counts, all_context_ids] = get_allcounts(eptrials, clusters, task_context_ids, window_duration, pre_time, iti_time, post_time, min_conf)

        %preallocate
        all_counts = cell(length(task_context_ids), 2);
        all_context_ids = cell(length(task_context_ids), 2);
        %befores then afters
        for firsthalf_2ndhalf = 1:2

            %iterate through contexts
            for ctxt_counter = 1:length(task_context_ids)
                current_context = task_context_ids(ctxt_counter);

                %prepost or iti time window
                if firsthalf_2ndhalf == 1 && ctxt_counter == 1
                    kt_tw = pre_time;
                elseif firsthalf_2ndhalf == 2 && ctxt_counter == 4
                    kt_tw = post_time;
                else
                    kt_tw = iti_time;
                end
                
                %key time
                if firsthalf_2ndhalf ==1
                    keytime = min(eptrials(eptrials(:,5)==ctxt_counter, 1));
                    keytime_start = keytime - kt_tw;
                    keytime_end = keytime+(6*60);
                elseif firsthalf_2ndhalf ==2
                    keytime = max(eptrials(eptrials(:,5)==ctxt_counter, 1));
                    keytime_start = keytime -(6*60);  
                    keytime_end = keytime + kt_tw;
                end
                keytime_idx = eptrials(:,1) > keytime_start  &  eptrials(:,1) < keytime_end;
                
                %calculate spike counts
                [spike_counts, tw_contexts] = slide_rate_window(eptrials(keytime_idx,:), clusters(clusters(:,2)>=min_conf,1), window_duration);

                %load into context-appropriate cell
                all_counts{ctxt_counter, firsthalf_2ndhalf} = [all_counts{ctxt_counter, firsthalf_2ndhalf} spike_counts];
                all_context_ids{ctxt_counter, firsthalf_2ndhalf} = [all_context_ids{ctxt_counter, firsthalf_2ndhalf} tw_contexts];
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


%pca figure function
%
    function comb_fig(rm, ctx_idx, first_second_idx, visit_mtx)
        
        %color
        colors = [0.6350    0.0780    0.1840;...
                  0.8500    0.3250    0.0980;...
                  0.9290    0.6940    0.1250;...
                  0.4660    0.6740    0.1880;...
                  0.1961    0.1961    0.1961;...
                  0.4510    0.4510    0.4510;...
                  0.6039    0.6039    0.6039;];
        
        %principle components
        [rm_pca, ~, ~, ~, explained_pca] = pca(rm');
        explained_pca = sum(explained_pca(1:3));

        rm_pca_vis = rm_pca(visit_mtx(:,3)==2, :);
        ctx_idx_vis = ctx_idx(visit_mtx(:,3)==2);
        rm_pca_iti = rm_pca(visit_mtx(:,3)~=2, :);
        ctx_idx_iti = visit_mtx(visit_mtx(:,3)~=2, 2);

        
        %plot
        figure; hold on
        count = 0; %visit tws
        for ctxti = unique(ctx_idx)'
            count = count+1;
            plot(rm_pca_vis(ctx_idx_vis==ctxti,1), rm_pca_vis(ctx_idx_vis==ctxti,2), '.', 'MarkerSize', 8, 'Color', colors(count,:));
        end
        
        count = 0; %iti tws
        iti_count = 0;
        for ctxti = unique(ctx_idx_iti)'
            count = count+1;
            if ismember(count, 2:4)
                iti_count = iti_count+1;
                plot(rm_pca_iti(ctx_idx_iti==ctxti,1), rm_pca_iti(ctx_idx_iti==ctxti,2), '.', 'MarkerSize', 8, 'Color', colors(iti_count+4,:));
            end
        end
        

        %{
        count = 0; %means
        for ctxti = unique(ctx_idx)'
            count = count+1
            h{count} = plot(mean(rm_pca_vis(ctx_idx_vis==ctxti,1)), mean(rm_pca_vis(ctx_idx_vis==ctxti,2)), '.', 'MarkerSize', 90, 'Color', colors(count,:)); 
            plot(mean(rm_pca_iti(ctx_idx_iti==ctxti,1)), mean(rm_pca_iti(ctx_idx_iti==ctxti,2)), '.', 'MarkerSize', 90, 'Color', colors(count,:)); 
        end
        %}
        %{
        xlabel('pc1')
        ylabel('pc2')
        zlabel('pc3')
        title('zscored spike counts')
        legend show
        legend([h{1}, h{2}, h{3}, h{4}], 'Visit1', 'Visit2', 'Visit3', 'Visit4')
%}
        axis square
        set(gca,'TickLength',[0, 0]); box off


    end


%redline function
%
    function [visit_mtx, ent_ext] = redlines(visit_mtx, varargin)
    %plot red lines using changes in context ids (visit_mtx(:,2))  
        
        if nargin == 2
            plot_on = varargin{1};
        else
            plot_on = 1;
        end
    
        %add two columns (ctxt or iti, visit number)
        visit_mtx = [visit_mtx nan(size(visit_mtx,1),2)];
        
        previous_ctx = visit_mtx(1,2); 
        iti_or_vis = 1;
        current_vis = 0;
        
        %transition times
        ent_ext = nan(4,2);

        for itw = 1:size(visit_mtx,1) %iter through ctx ids

            current_ctx = visit_mtx(itw,2);
            
            if current_ctx == previous_ctx %if same

                visit_mtx(itw,3) = iti_or_vis;
                
                if iti_or_vis ==2
                    visit_mtx(itw,4) = current_vis;
                end
                
                continue

            elseif iti_or_vis == 1

                if ~ismember(current_ctx, 9:13) && ~isnan(current_ctx) %not iti

                    if plot_on == 1
                        hold on 
                        plot([visit_mtx(itw,1)-.5 visit_mtx(itw,1)-.5], [0 10], 'r-', 'linewidth', 1)
                    end
                    
                    
                    iti_or_vis = 2;
                    current_vis = current_vis + 1;
                    ent_ext(current_vis, 1) = visit_mtx(itw,1);
                    visit_mtx(itw,3:4) = [iti_or_vis current_vis];
                    previous_ctx = visit_mtx(itw,2);
                else
                    visit_mtx(itw,3) = iti_or_vis;
                    previous_ctx = visit_mtx(itw,2);
                end

            elseif iti_or_vis == 2

                if ~ismember(current_ctx, 1:8) %not visit

                    if plot_on == 1
                        hold on 
                        plot([visit_mtx(itw,1)-.5 visit_mtx(itw,1)-.5], [0 10], 'r-', 'linewidth', 1)
                    end
                    
                    iti_or_vis = 1;
                    ent_ext(current_vis, 2) = visit_mtx(itw,1);
                    visit_mtx(itw,3) = iti_or_vis;
                    previous_ctx = visit_mtx(itw,2);
                else
                    visit_mtx(itw,3:4) = [iti_or_vis current_vis];
                    previous_ctx = visit_mtx(itw,2);
                end
            end
        end
        
        
    end


%fitline function
%{
    function [R, p] = fit_line(X,Y)
        %fits a line to a scatterplot of X,Y
        poly=polyfit(X(~isnan(X) & ~isnan(Y)),Y(~isnan(X) & ~isnan(Y)), 1);        
        %fit_x = (min(X)-(range(X)/12)):(range(X)/100):(max(X)+(range(X)/12));
        X_hold = X(~isnan(X) & ~isnan(Y));
        fit_x = (min(X_hold)-(range(X_hold)/10)):(range(X_hold)/100):(max(X_hold)+(range(X_hold)/10));
        fit_y=polyval(poly, fit_x);
        
        hold on; 
        plot(fit_x, fit_y, 'LineWidth', 2, 'Color', [.5 .5 .5])

        [R, p] = corr(X(~isnan(X) & ~isnan(Y)),Y(~isnan(X) & ~isnan(Y)));

    end
%}





end