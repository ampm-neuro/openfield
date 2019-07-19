function shuffles_out = shuffle_corr(mtx, num_shuffles)
%shuffles all_rates to find out how correllated the context representations
%were by chance. Chance determined by shuffling individual cells between
%contexts.


mtx_size = size(mtx);

shuffles_out = nan(num_shuffles, 4);



for i = 1:num_shuffles
    
    %shuffle context id of each cell rate
    shuffle_idx_1 = nan(mtx_size);
    shuffle_idx_2 = nan(mtx_size);
    for s = 1:mtx_size(1)
        shuffle_idx_1(s, :) = randperm(mtx_size(2));
        shuffle_idx_2(s, :) = randperm(mtx_size(2));
    end
    mtx_shuf_1 = mtx(shuffle_idx_1);
    mtx_shuf_2 = mtx(shuffle_idx_2);

    %corm
    cm_whole = corm(mtx_shuf_1, mtx_shuf_2);
    shuffles_out(i,:) = [cm_whole(1,7) mean([cm_whole(1,2); cm_whole(2,7)]) cm_whole(3,4) cm_whole(5,6)];%bookend, b&w, obj, arrang
end