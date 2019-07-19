function rate_mtx = collapse_time_individual(rate_mtx, context_id, visit_order)
%takes ouput of all_countwindow (all_counts) and collapses across the
%time dimension 

       ic = 10;
       %figure; imagesc(rate_mtx); title A_rm
       
       %figure; plot(smooth(zscore_mtx(rate_mtx(:,ic)),10)); title A_rm
      
       %number every time window
       tw_ct = 1:size(rate_mtx,1);
       
       %sort rate matrix based on visit order
       %
       
       %rename context visits and itis to be easily sortable
       old_assign = [1 2 3 4];
       new_assign = [11 13 15 17];
       visit_order_new = visit_order;
       count = 0;
       for v_old = old_assign
           count = count+1;
           visit_order_new(visit_order_new==v_old) = new_assign(count);
       end
        %also itis
        context_id = repmat(context_id, 1, size(visit_order,2));
        old_assign = [10 11 12];
        new_assign = [12 14 16];
        count = 0;
        for v_old = old_assign
           count = count+1;
           visit_order_new(context_id==v_old) = new_assign(count);
        end
       
       %figure; fit_line(1:length(zscore_mtx(rate_mtx(:,ic))), zscore_mtx(rate_mtx(:,ic))); title a_rm
        
       %sort
       [~, sorted_von_idx] = sort(visit_order_new);
       unsort_idx_prep = repmat((1:size(sorted_von_idx,1))', 1, size(sorted_von_idx,2));
       unsort_idx = nan(size(unsort_idx_prep));
       sorttest = unsort_idx_prep;
    
       sorted_rm = rate_mtx;
       for col = 1:size(sorted_von_idx,2)
           sorted_rm(:,col) = rate_mtx(sorted_von_idx(:, col),col);
           sorttest(:,col) = unsort_idx_prep(sorted_von_idx(:, col),col);
           unsort_idx(sorted_von_idx(:, col), col) = unsort_idx_prep(:, col);
       end
       
       %figure; fit_line(1:length(zscore_mtx(sorted_rm(:,ic))), zscore_mtx(sorted_rm(:,ic))); title b_rm
              
       %remove effect of time on firing rate
       for icell = 1:size(sorted_rm,2)
           p = polyfit(tw_ct', sorted_rm(:,icell), 1);
           sorted_rm(:,icell) = sorted_rm(:,icell) - p(1).*(tw_ct');
       end
       
       %smooth
       for i = 1:size(sorted_rm,2)
           sorted_rm(:,i) = smooth(sorted_rm(:,i),5);
       end
       
       %figure; fit_line(1:length(zscore_mtx(sorted_rm(:,ic))), zscore_mtx(sorted_rm(:,ic))); title c_rm
              
       %sort back
       for col = 1:size(unsort_idx,2)
           rate_mtx(:,col) = sorted_rm(unsort_idx(:, col),col);
           sorttest(:,col) = sorttest(unsort_idx(:, col),col);
       end
         
       %figure; plot(smooth(zscore_mtx(rate_mtx(:,ic)),10)); title B_rm
       %figure; fit_line(1:length(zscore_mtx(rate_mtx(:,ic))), zscore_mtx(rate_mtx(:,ic))); title d_rm
       
 end