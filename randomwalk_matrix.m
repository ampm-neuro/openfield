%
%preallocate random walk rate matrix
rm_rw = nan(size(rm));

%for all cells
for ic = 1:size(rm, 2)

    %typical rate changes between timewindows
    typical_rc = rm(2:end, ic)-rm(1:end-1, ic);
    
    %remove any positive or negative bias
    %TEST WITH AND WITHOUT
    posneg = rand(size(typical_rc))<=.5;
    typical_rc = typical_rc.*posneg;
    
    %randomize order of rate changes
    typical_rc = typical_rc(randperm(length(typical_rc)));
    
    %starting point same for shuf and orig
    rm_rw(1,ic) = rm(1,ic); 

    %cumulative rate changes
    rm_rw(:,ic) = cumsum([rm(1,ic) typical_rc']);

end

%zscore
rm_rw_z = zscore_mtx(rm_rw);
%}

%Class and plot
%{
[p_correct_v1_shufrw, assignments_v1_shufrw, distances_v1_shufrw] = ...
    md_classify(rm_rw_z(visit_idx(:,4)==1, :), visit_idx(visit_idx(:,4)==1,5));
figure 
imagesc(histcounts2(visit_idx(visit_idx(:,4)==4,5), assignments_v1_shufrw, 24)) 
%}
%
[p_correct_shufrw, assignments_shufrw, distances_shufrw] = ...
    md_classify(rm_rw_z, visit_idx(:,5));
figure 
imagesc(histcounts2(visit_idx(:,5), assignments_v1_shufrw)) 
%}

caxis([0 40]);axis square; colorbar