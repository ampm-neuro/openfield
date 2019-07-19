%load('pop_time_prep')

tic
num_shufs = 100;
pshufs = nan(num_shufs, 4);

%shuffle
for shuf = 1:num_shufs

    %shuffle each cell's rates
    rm_z_shuf = nan(size(rm_z));
    for icol = 1:size(rm_z,2) 
        hold = rm_z(:,icol); 
        rm_z_shuf(:,icol) = hold(randperm(length(hold))); 
    end

    %caculate each visit's classification accuracy
    for vis = 1:4
        pshufs(shuf, vis) = md_classify(rm_z_shuf(visit_idx(:,4)==vis,:), visit_idx(visit_idx(:,4)==vis,5));
    end
    
end
toc