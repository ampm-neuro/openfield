%pull out populations
load('C:\Users\ampm1\Desktop\oldmatlab\PCia\pop_ctx_9_27_simple.mat')
load('C:\Users\ampm1\Desktop\oldmatlab\PCia\250ms_with_hd.mat')
%fraction of interest
foi = 1/2;

all_distances = [];
% (remove columns)
posterior_for_vid = [];
bw_pop_cell_idx_hold = bw_pop_cell_idx;
spk_cts = all_counts_bw;
spk_cts = spk_cts(:, bw_pop_cell_idx_hold > 0);
posx = all_hd_bw(:,bw_pop_cell_idx_hold > 0,1);
bw_pop_cell_idx_hold = bw_pop_cell_idx_hold(bw_pop_cell_idx_hold > 0);

%cut timewindows from irrelevant contexts
% (remove rows)
rctxts = unique(context_bw(context_bw<9))'; %1:4
rctxts_idx = ismember(context_bw, rctxts);
spk_cts = spk_cts(rctxts_idx, :);
posx = posx(rctxts_idx, :);
ctxt_idx = context_bw(rctxts_idx, :);


%timespace bin index
%

%find spatial bin classifications
%0 and 360 are hard coded min and max for both x and y (see trials_pcia)
%bins = 6;
bins = 7; bins = bins^2;
bin_edges = linspace(0, 360, bins+1);

%error prevention
posx(posx<0) = 0;
posx(posx>360) = 360;


%context of interest
CoI = 3;

%for each context of interest
for ictxt = CoI

    %for each cell (shorthand for session -> could be improved)
    for iclust = 1:size(spk_cts,2)

        %find bin indices for each x and... 
        posx_hold = posx(ctxt_idx==ictxt, iclust);
        [~, ~, bin_idx_x] = histcounts(posx_hold, bin_edges);
        bin_idx_x(isnan(posx_hold)) = nan;

        %exclude nans
        rng_idx = bin_idx_x>0 & bin_idx_x<=bins;

        %load space index by context and cluster
        timespace_idx(ctxt_idx==ictxt, iclust) = bin_idx_x(rng_idx); 

    end
end

%preallocate context_matchup X population  
out_class_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
out_true_wc = cell(1/foi, length(unique(bw_pop_cell_idx_hold)));
accuracy_wc = nan(1/foi, length(unique(bw_pop_cell_idx_hold)));

%bayesian decode each pop session
for ipop = unique(bw_pop_cell_idx_hold)'
    
    %context spike counts
    spike_counts = spk_cts(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
    %context classifications
    classifications = timespace_idx(ctxt_idx == CoI, bw_pop_cell_idx_hold == ipop);
    classifications = classifications(:,1); %columns are redundant
    
    %fill cell for iteratations below
    %fraction of interest
    train_test = cell(1/foi, 2);
    for i = 1:1/foi
        foi_indx = floor(length(classifications)/(1/foi)*(i-1)+1):floor(length(classifications)/(1/foi)*i);
        train_test{i,1} = spike_counts(foi_indx, :);
        train_test{i,2} = classifications(foi_indx);
    end

    
    %for each half
    for iwc = 1:1/foi

        %training sample (context)
        training_wc = cat(1,train_test{setdiff(1:1/foi, iwc), 1});
       
        %true training class
        group_train_wc = cat(1,train_test{setdiff(1:1/foi, iwc), 2});

        %track excluded pixles
        included_pixels = 1:bins;
        
        %remove training pixles that weren't visited for at least 5sec.
        %
        [a, b] = histc(group_train_wc, unique(group_train_wc));
        deletion_idx = ismember(b, find(a<(5*4)));
        group_train_wc(deletion_idx) = [];
        training_wc(deletion_idx, :) = [];
        included_pixels = intersect(included_pixels, unique(group_train_wc));
        %}
        
        %testing sample (context)
        sample_wc = train_test{iwc, 1};
        
        
        
        %true test class
        group_test_wc = train_test{iwc, 2}; 

        %remove test pixles that aren't in train set
        %
        deletion_idx_test = ...
            ~ismember(group_test_wc, group_train_wc);
        sample_wc(deletion_idx_test, :) = [];
        group_test_wc(deletion_idx_test) = [];
        %}
        
        %reliability plots
        %{
        if ipop==5
           
           f_x = nan(length(unique(group_train_wc)), size(training_wc,2));
           for pixel = 1:size(f_x,1)
                f_x(pixel, :) = mean(training_wc(group_train_wc==pixel, :))/0.25;
           end

           if iwc == 1
            fig1 = figure; plot(f_x(:,1)); title('cell 1'); set(gca,'TickLength',[0, 0]);
            fig2 = figure; plot(f_x(:,2)); title('cell 2'); set(gca,'TickLength',[0, 0]);
            fig3 = figure; plot(f_x(:,3)); title('cell 3'); set(gca,'TickLength',[0, 0]);
            fig4 = figure; plot(f_x(:,4)); title('cell 4'); set(gca,'TickLength',[0, 0]);
            fig5 = figure; plot(f_x(:,5)); title('cell 5'); set(gca,'TickLength',[0, 0]);
            fig6 = figure; plot(f_x(:,6)); title('cell 6'); set(gca,'TickLength',[0, 0]);
            fig7 = figure; plot(f_x(:,7)); title('cell 7'); set(gca,'TickLength',[0, 0]);
            fig8 = figure; plot(f_x(:,8)); title('cell 8'); set(gca,'TickLength',[0, 0]);
            fig9 = figure; plot(f_x(:,9)); title('cell 9'); set(gca,'TickLength',[0, 0]);
            fig10 = figure; plot(f_x(:,10)); title('cell 10'); set(gca,'TickLength',[0, 0]);
            fig11 = figure; plot(f_x(:,11)); title('cell 11'); set(gca,'TickLength',[0, 0]);


           elseif iwc == 1/foi
               figure(fig1); hold on; plot(f_x(:,1)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                                var_name = 'cell1'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')
               figure(fig2); hold on; plot(f_x(:,2)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell2'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig3); hold on; plot(f_x(:,3)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell3'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

              figure(fig4); hold on; plot(f_x(:,4)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                             var_name = 'cell4'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig5); hold on; plot(f_x(:,5)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell5'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig6); hold on; plot(f_x(:,6)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell6'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

              figure(fig7); hold on; plot(f_x(:,7)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                             var_name = 'cell7'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig8); hold on; plot(f_x(:,8)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell8'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig9); hold on; plot(f_x(:,9)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell9'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

              figure(fig10); hold on; plot(f_x(:,10)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                             var_name = 'cell10'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

               figure(fig11); hold on; plot(f_x(:,11)); hold off; xticks([0:5:bins+1]);xticklabels(round(xticks.*(360/bins)))
                              var_name = 'cell11'; print(['C:\Users\ampm1\Desktop\Maze_Revisions\HD\' var_name], '-dpdf', '-painters', '-bestfit')

           else
               figure(fig1); hold on; plot(f_x(:,1)); hold off; 
               figure(fig2); hold on; plot(f_x(:,2)); hold off; 
               figure(fig3); hold on; plot(f_x(:,3)); hold off; 
              figure(fig4); hold on; plot(f_x(:,4)); hold off; 
               figure(fig5); hold on; plot(f_x(:,5)); hold off; 
               figure(fig6); hold on; plot(f_x(:,6)); hold off;
              figure(fig7); hold on; plot(f_x(:,7)); hold off; 
               figure(fig8); hold on; plot(f_x(:,8)); hold off; 
               figure(fig9); hold on; plot(f_x(:,9)); hold off;
              figure(fig10); hold on; plot(f_x(:,10)); hold off; 
               figure(fig11); hold on; plot(f_x(:,11)); hold off;
               
           end
           
        end
        %}
        
        
        
        %decode
        [class_wc, posterior_wc] = bayesian_decode(sample_wc, training_wc, group_train_wc, .25);

        %posterior_hold(:, included_pixels) = posterior_wc;
        posterior_wc = posterior_wc./sum(posterior_wc(:));
        
        %load
        out_class_wc{iwc, ipop} =  class_wc;
        out_true_wc{iwc, ipop} =  group_test_wc;
        accuracy_wc(iwc, ipop) = sum(class_wc==group_test_wc)/length(group_test_wc);
        
    end
    
    %distances (in degrees) between decoded positions and true positions
    all_class = [cell2mat(out_class_wc(1,ipop)); cell2mat(out_class_wc(2,ipop))];
    all_true = [cell2mat(out_true_wc(1,ipop)); cell2mat(out_true_wc(2,ipop))];

        for i = 1:length(all_true)
            all_distances = [all_distances; circ_distance(all_class(i), all_true(i), [1 bins])];
        end
    
    accuracy_distances(ipop) = sum(all_distances <= (100/bins)*sqrt(2) )/length(all_distances);
   
end


%plot
%
figure; hold on
cell_e{1} = mean(accuracy_wc)';
errorbar_plot(cell_e)
ylim([0 inf])
plot([0.5 2.5], [1/bins 1/bins], 'k--')
title absolute
%}
%
figure; hold on
cell_e{1} = accuracy_distances;
errorbar_plot(cell_e)
ylim([0 inf])
title distances
%}





