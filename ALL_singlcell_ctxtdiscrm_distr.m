function [rate_cell_out, space_cell_out, all_ctxt_rate_b, all_ctxt_rate_w, bw_zs, bw_ts, ai_ts] = ALL_singlcell_ctxtdiscrm_distr
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

all_ctxt_rate_b = [];
all_ctxt_rate_w = [];
bw_zs = [];
bw_ts = [];
ai_ts = [];

tw_duration = 30;%s

discrim_distr_bw_rate = [];
discrim_distr_bw_rate_within = [];
discrim_distr_a_rate = [];
discrim_distr_i_rate = [];
discrim_distr_ia_rate = [];


distr_bw_within_space = [];
distr_bw_between_space = [];
discrim_distr_bw_space = [];
%discrim_distr_a_space = [];
%discrim_distr_i_space = [];
%discrim_distr_ia_space = [];

%names%get all the things in neurodata folder...
file_list_subjects = dir('neurodata\');

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
    rat = file_names_subjects{:}(subject,1).name
    
    %task folders
    task_folders = {'black_and_white' 'object_arrangement'};
    
    for task = 1:2
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\', num2str(file_names_subjects{:}(subject,1).name),'\', task_folders{task}, '\*.mat'));


        %number of folders
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name
            %day = session

            %load session
            eptrials = [];
            context_ids = [];
            origin_file = [];
            load(strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file)));
            strcat('C:\Users\ampm1\Desktop\oldmatlab\PCia\neurodata\',num2str(rat),'\', task_folders{task},'\' ,num2str(day_file))
            
            %if ~isequal(context_ids, [3 1 2 4]')
            %    continue
            %end
            
            if isempty(clusters)
                continue
            else
                %clusters = clusters(clusters(:,2)>2, 1);
                clusters = clusters(:, 1);
                if isempty(clusters)
                    continue
                end
            end

            %calculate mean rate from each context visit
            [full_visit, fh_visit, sh_visit] = ctxt_rates(eptrials, clusters, sort(context_ids));
            %normalize by sum
            full_visit = full_visit./sum(full_visit,2);
            fh_visit = fh_visit./sum([fh_visit sh_visit],2);
            sh_visit = sh_visit./sum([fh_visit sh_visit],2);
            
            %mean rate in each pixle from each context visit
            [vectors_first, vectors_second, vectors_all] = pop_pixlecorr(eptrials, sort(context_ids), clusters, 40);
            

            %If BW session
            if isequal(sort(context_ids), [1 2 3 4]')

                %FULL SESSIONS
                %
                    %rate
                    %{
                    vis_id = [];
                    rates = [];
                    for ctxt = 1:4
                        scs{ctxt} = slide_rate_window(eptrials(eptrials(:,6)==ctxt, :), clusters(:,1), tw_duration);
                        scs{ctxt} = scs{ctxt}./tw_duration;
                        %visidx{ctxt} = ones(size(scs{ctxt})).*mode(eptrials(eptrials(:,6)==ctxt, 5)); %keep track of visit number
                    
                        vis_id = [vis_id ; ones(size(scs{ctxt},1),1).*mode(eptrials(eptrials(:,6)==ctxt, 5))];
                        rates = [rates ; scs{ctxt}./tw_duration];
                    end 
                    
                    %figure; 
                    %subplot(1,2,1)
                    %plot(scs{1}(:,1))
                   
                    
                    %remove effect of time on firing rate
                    rt_out = removetime(vis_id, rates);
                    %}
                    [rates_all, ctx_ids] = slide_rate_window(eptrials, clusters(:,1), tw_duration);
                    rates_all = zscore_mtx(rates_all);
                    
                    
                    %put back into context id cells
                    count = 1;
                    for ctxt = 1:4
                        %scs{ctxt} = rt_out(count:count+length(scs{ctxt})-1, :);
                        %count = count+length(scs{ctxt});
                        
                        scs{ctxt} = rates_all(ctx_ids == ctxt, :);
                    end
                    
                    %subplot(1,2,2)
                    %plot(scs{1}(:,1))
                    %figure; imagesc(rt_out)

                    %
                    between_full_rate = mean([...
                        (mean(scs{1}) - mean(scs{3}))    ./ mean([std(scs{1}); std(scs{3})]);  ...
                        (mean(scs{2}) - mean(scs{4}))    ./ mean([std(scs{2}); std(scs{4})])   ...
                        ])';
                    
                    within_full_rate = mean([...
                        (mean(scs{1}) - mean(scs{2}))    ./ mean([std(scs{1}); std(scs{2})]);  ...
                        (mean(scs{3}) - mean(scs{4}))    ./ mean([std(scs{3}); std(scs{4})])   ...
                        ])';
                    %}
                    
                    
                    %(mean(black_ctx_rates) - mean(white_ctx_rates))';
                    %between_full_rate_std = mean([std(black_ctx_rates); std(white_ctx_rates)])';
                    %between_full_rate_std = mean([std(scs{1}); std(scs{2}); std(scs{3}); std(scs{4})])';
                    for ic = 1:size(clusters,1)
                        if ~isequal(length(scs{1}(:,ic)), length(scs{2}(:,ic)), length(scs{3}(:,ic)), length(scs{4}(:,ic)))
                            minlen = min([length(scs{1}(:,ic)), length(scs{2}(:,ic)), length(scs{3}(:,ic)), length(scs{4}(:,ic))]);
                            scs1 = scs{1}(1:minlen,ic);
                            scs2 = scs{2}(1:minlen,ic);
                            scs3 = scs{3}(1:minlen,ic);
                            scs4 = scs{4}(1:minlen,ic);
                        else
                            scs1 = scs{1}(:,ic);
                            scs2 = scs{2}(:,ic);
                            scs3 = scs{3}(:,ic);
                            scs4 = scs{4}(:,ic);
                        end
                        
                        Y = [scs1; scs2; scs3; scs4];
                        S = [(1:length(scs1))' ;  (1:length(scs2))' ; (1:length(scs3))'; (1:length(scs4))' ];
                        F1 = [ones(size(scs1)); ones(size(scs2)).*2; ones(size(scs3)); ones(size(scs4)).*2];
                        F2 = [ones(size([scs1; scs2])); ones(size([scs3; scs4])).*2];

                        stats = rm_anova2(Y,S,F1,F2,{'vis_12', 'ctx_bw'});

                        between_full_rate_anova = stats{3,6};
                        within_full_rate_anova = stats{2,6};

                        all_ctxt_rate_b = [all_ctxt_rate_b; between_full_rate_anova];
                        all_ctxt_rate_w = [all_ctxt_rate_w; within_full_rate_anova];
                        
                        
                        
                        bz = (mean([scs1; scs2]) - mean([scs3; scs4]))/std(Y);
                        wz = (mean([scs1; scs3]) - mean([scs2; scs4]))/std(Y);
                        bw_zs = [bw_zs; [bz wz]];
                        
                        [~,~,~, bt] = ttest([scs1; scs2],[scs3; scs4]);
                        [~,~,~, wt] = ttest([scs1; scs3],[scs2; scs4]);
                        bw_ts = [bw_ts; [bt.tstat wt.tstat]];
                        
                        %PLOT EXAMPLE CELLS
                        %{
                        if bt.tstat < -5 && abs(wt.tstat) < 2
                            heatmaps_trial(eptrials, clusters(ic), [3 1 4 2]);
                        end
                        %}
                        
                        
                    end
                    
                    
                    %space correlation
                    within_full_corr_space = mean([...
                                  diag(corr(vectors_all(:,:,1),vectors_all(:,:,2)))...
                                  diag(corr(vectors_all(:,:,3),vectors_all(:,:,4)))...
                                                           ],2);
                    between_full_corr_space = mean([...
                                  diag(corr( vectors_all(:,:,1),vectors_all(:,:,3) ))...
                                  diag(corr( vectors_all(:,:,1),vectors_all(:,:,4) ))...
                                  diag(corr( vectors_all(:,:,2),vectors_all(:,:,3) ))...
                                  diag(corr( vectors_all(:,:,2),vectors_all(:,:,4) ))...
                                                           ],2);
                
                
                %HALF SESSIONS
                %
                %between_full_rate_meandiff = abs(mean(full_visit(:, [1 2]),2) - mean(full_visit(:, [3 4]),2));
                %between_full_rate_std = mean(std(full_visit(:, [1 2]),1,2) - std(full_visit(:, [3 4]),1,2));
                
                
                %OUTPUT
                %
                
                %discrim index         
                discrim_distr_sesh_rate = between_full_rate;%./between_full_rate_std;
                discrim_distr_sesh_rate_within = within_full_rate;
                
                discrim_distr_sesh_space = within_full_corr_space - between_full_corr_space;
                
                
                
                %{
                zscore_rate_diff = 3;
                if any(abs(discrim_distr_sesh_rate) > zscore_rate_diff)
                    
                    strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file))
                    
                    for clust = clusters(abs(discrim_distr_sesh_rate) > zscore_rate_diff,1)'
                    
                        heatmaps_trial(eptrials, clust, context_ids);
                    
                        figure; hold on
                        title(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)))
                        for ctxt_vis = 1:4
                            
                            c_ctxt = context_ids(ctxt_vis);
                            
                            scs{ctxt_vis} = slide_rate_window(eptrials(eptrials(:,6)==c_ctxt, :), clust, tw_duration);
                            scs{ctxt_vis} = scs{ctxt_vis}./tw_duration;
                            plot(repmat(ctxt_vis, size(scs{ctxt_vis})), scs{ctxt_vis}, 'o', 'markersize', 9, 'color', [.8 .8 .8])
                        end
                        errorbar([mean(scs{1}) mean(scs{2}) mean(scs{3}) mean(scs{4})],...
                            [std(scs{1})/sqrt(length(scs{1})) std(scs{2})/sqrt(length(scs{2})) std(scs{3})/sqrt(length(scs{3})) std(scs{4})/sqrt(length(scs{4}))], 'k')
                        xlim([.5 4.5]); set(gca,'TickLength',[0, 0]); box off
                    end
                end
                %}
                %{
                corr_within = 0.2;
                corr_between = 0;
                %corr_diff = 0.3;
                if any(within_full_corr_space > corr_within) && any(between_full_corr_space < corr_between)
                    
                    strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file))
                    
                    
                    for clust = clusters(within_full_corr_space > corr_within & between_full_corr_space < corr_between,1)'
                        
                        heatmaps_trial(eptrials, clust, context_ids);

                        figure; hold on
                        title(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)))
                        for ctxt_vis = 1:4

                            c_ctxt = context_ids(ctxt_vis);

                            scs{ctxt_vis} = slide_rate_window(eptrials(eptrials(:,6)==c_ctxt, :), clust, tw_duration);
                            scs{ctxt_vis} = scs{ctxt_vis}./tw_duration;
                            plot(repmat(ctxt_vis, size(scs{ctxt_vis})), scs{ctxt_vis}, 'o', 'markersize', 9, 'color', [.8 .8 .8])
                        end
                        errorbar([mean(scs{1}) mean(scs{2}) mean(scs{3}) mean(scs{4})],...
                            [std(scs{context_ids(1)})/sqrt(length(scs{1})) std(scs{2})/sqrt(length(scs{2})) std(scs{3})/sqrt(length(scs{3})) std(scs{4})/sqrt(length(scs{4}))], 'k')
                        xlim([.5 4.5]); set(gca,'TickLength',[0, 0]); box off
                    end

                end
                %}
                
                
                %load
                discrim_distr_bw_rate = [discrim_distr_bw_rate; discrim_distr_sesh_rate];
                discrim_distr_bw_rate_within = [discrim_distr_bw_rate_within; discrim_distr_sesh_rate_within];
                distr_bw_within_space = [distr_bw_within_space; within_full_corr_space];
                distr_bw_between_space = [distr_bw_between_space; between_full_corr_space];
                discrim_distr_bw_space = [discrim_distr_bw_space; discrim_distr_sesh_space];
                
                
                
                
            %if IA session
            %
            elseif isequal(sort(context_ids), [5 6 7 8]')
                
                ia_ctxts = 5:8;
                
                %rate
                vis_id = [];
                rates = [];
                
                [rates_all, ctx_ids] = slide_rate_window(eptrials, clusters(:,1), tw_duration);
                rates_all = zscore_mtx(rates_all);
                %{
                for ctxt = 1:4
                    scs{ctxt} = slide_rate_window(eptrials(eptrials(:,6)==ia_ctxts(ctxt), :), clusters(:,1), tw_duration);
                    scs{ctxt} = scs{ctxt}./tw_duration;

                    vis_id = [vis_id ; ones(size(scs{ctxt},1),1).*mode(eptrials(eptrials(:,6)==ia_ctxts(ctxt), 5))];
                    rates = [rates ; scs{ctxt}./tw_duration];
                end 
                %}

                %remove effect of time on firing rate
                %
                %rt_out = removetime(vis_id, rates);
                %rt_out = removetime(1:length(rates_all), rates_all);


                %put back into context id cells
                count = 1;
                for ctxt = 1:4
                    %scs{ctxt} = rt_out(count:count+length(scs{ctxt})-1, :);
                    %count = count+length(scs{ctxt});
                    
                    scs{ctxt} = rates_all(ctx_ids == ia_ctxts(ctxt), :);
                    
                    
                end
                %}


                for ic = 1:length(clusters)
                    if ~isequal(length(scs{1}(:,ic)), length(scs{2}(:,ic)), length(scs{3}(:,ic)), length(scs{4}(:,ic)))
                        minlen = min([length(scs{1}(:,ic)), length(scs{2}(:,ic)), length(scs{3}(:,ic)), length(scs{4}(:,ic))]);
                        scs1 = scs{1}(1:minlen,ic);
                        scs2 = scs{2}(1:minlen,ic);
                        scs3 = scs{3}(1:minlen,ic);
                        scs4 = scs{4}(1:minlen,ic);
                    else
                        scs1 = scs{1}(:,ic);
                        scs2 = scs{2}(:,ic);
                        scs3 = scs{3}(:,ic);
                        scs4 = scs{4}(:,ic);
                    end

                    [~,~,~, bt] = ttest(scs1, scs2);
                    [~,~,~, wt] = ttest(scs3, scs4);
                    ai_ts = [ai_ts; [bt.tstat wt.tstat]]; %arrangement, object

                    %PLOT EXAMPLE CELLS
                    %{
                    if bt.tstat < -5 && abs(wt.tstat) < 2
                        heatmaps_trial(eptrials, clusters(ic), [3 1 4 2]);
                    end
                    %}


                end
                    
                
                %rate
                %{
                for ctxt_idx = 1:4

                    scs{ctxt_idx} = slide_rate_window(eptrials(eptrials(:,6)==ia_ctxts(ctxt_idx), :), clusters(:,1), tw_duration);
                    scs{ctxt_idx} = scs{ctxt_idx}./tw_duration;
                    
                end 
                
                %discrim index                    
                discrim_distr_sesh_a = (mean(scs{1}) - mean(scs{2})) ./ mean([std(scs{1}); std(scs{2})]);
                discrim_distr_sesh_i = (mean(scs{3}) - mean(scs{4})) ./ mean([std(scs{3}); std(scs{4})]);
                
                %discrim_distr_sesh_ia = (between_diff_ia-within_diff_ia)./(within_diff_ia+between_diff_ia);

                %load
                discrim_distr_a_rate = [discrim_distr_a_rate; discrim_distr_sesh_a'];
                discrim_distr_i_rate = [discrim_distr_i_rate; discrim_distr_sesh_i'];
                %discrim_distr_ia_rate = [discrim_distr_ia_rate; discrim_distr_sesh_ia];
                %{
                zscore_diff = 2;
                if any(abs(discrim_distr_sesh_a) > zscore_diff)
                    
                    strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file))
                    
                    
                    for clust = clusters(abs(discrim_distr_sesh_a) > zscore_diff & abs(discrim_distr_sesh_a) <= 3, 1)'
                    
                        heatmaps_trial(eptrials, clust, context_ids);
                    
                        figure; hold on
                        title(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)))
                        for ctxt_vis = 1:4
                            
                            c_ctxt = context_ids(ctxt_vis);
                            
                            scs{ctxt_vis} = slide_rate_window(eptrials(eptrials(:,6)==c_ctxt, :), clust, tw_duration);
                            scs{ctxt_vis} = scs{ctxt_vis}./tw_duration;
                            plot(repmat(ctxt_vis, size(scs{ctxt_vis})), scs{ctxt_vis}, 'o', 'markersize', 9, 'color', [.8 .8 .8])
                        end
                        errorbar([mean(scs{1}) mean(scs{2}) mean(scs{3}) mean(scs{4})],...
                            [std(scs{1})/sqrt(length(scs{1})) std(scs{2})/sqrt(length(scs{2})) std(scs{3})/sqrt(length(scs{3})) std(scs{4})/sqrt(length(scs{4}))], 'k')
                        xlim([.5 4.5]); set(gca,'TickLength',[0, 0]); box off
                    end
                end
                %}
                %}
            else
                
                problem_ctxt_ids = context_ids
                continue
                %error('session type not detected')
            end
            %}
            
            
            
            
            %}

        end
    end
end

rate_cell_out{1} = discrim_distr_bw_rate;
rate_cell_out{2} = discrim_distr_a_rate;
rate_cell_out{3} = discrim_distr_i_rate;
rate_cell_out{4} = discrim_distr_ia_rate;
rate_cell_out{5} = discrim_distr_bw_rate_within;


space_cell_out{1} = distr_bw_within_space;
space_cell_out{2} = distr_bw_between_space;
space_cell_out{3} = discrim_distr_bw_space;
%space_cell_out{4} = discrim_distr_ia_rate;



 function rt_out = removetime(vis_id, rates)
                    
       %sort based on visit
       [sorted_rates, sorted_idx] = sortrows([vis_id rates], 1);
       sorted_rates = sorted_rates(:,2:end); %just cell rates
     
       %number every time window
       tw_ct = 1:size(sorted_rates,1);
       
       for icell = 1:size(sorted_rates,2)
           p = polyfit(tw_ct', sorted_rates(:,icell), 1);
           sorted_rates(:,icell) = sorted_rates(:,icell) - p(1).*(tw_ct');
       end
       [~, unsort_idx] = sort(sorted_idx);
       rt_out = sorted_rates(unsort_idx, :);
 end

end

