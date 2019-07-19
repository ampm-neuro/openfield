function shuffles_out = shuffle_recon(mtx1, mtx2, num_shuffles)
%shuffles for reconstruction error (first half, second half)


mtx_size = size(mtx1);

shuffles_out = nan(num_shuffles, 1);



for i = 1:num_shuffles
    
    %shuffle context id of each cell rate
    shuffle_idx_1 = nan(mtx_size);
    shuffle_idx_2 = nan(mtx_size);
    for s = 1:mtx_size(1)
        shuffle_idx_1(s, :) = randperm(mtx_size(2));
        shuffle_idx_2(s, :) = randperm(mtx_size(2));
    end
    mtx_shuf_1 = mtx1(shuffle_idx_1);
    mtx_shuf_2 = mtx2(shuffle_idx_2);

    %corm
    cm_halves = corm(mtx_shuf_1, mtx_shuf_2);
    
    %reconstruct
    reconst=nan(mtx_size(2),1);
    for i2 = 1:mtx_size(2)
        reconst(i2) = find(cm_halves(i2,:)==max(cm_halves(i2,:)));
    end

    %load recon error (how many reconstruction errors)
    errors = reconst-(1:mtx_size(2))';
    shuffles_out(i) = sum(errors~=0);
    
    
end

end