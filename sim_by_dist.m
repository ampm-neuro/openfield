function sim_by_dist(pvects)
%load pop-pxlecorr.mat
%computes the similarity (rval) between pixles as a function of distance

%number of pixels in pvects
n_px = size(pvects,1);
sides = sqrt(n_px);

%subsample pvects
%{
lo_hi = [0.3 0.8];
pvects_idx = reshape(1:n_px, sides, sides);
pvects_idx = pvects_idx(round(lo_hi(1)*sides : lo_hi(2)*sides), round(lo_hi(1)*sides : lo_hi(2)*sides));
pvects = pvects(pvects_idx(:), :,:);
n_px = size(pvects,1);
sides = sqrt(n_px);
%}

%compute similarity between all pixels
corr_mtx = nan(n_px^2, 4);
for ctxt = 1:(size(pvects,3)-1)
   corr_mtx(:,ctxt) = simvect(pvects, ctxt, n_px); 
end
[cm, pixel_ids] = simvect(pvects, size(pvects,3), n_px); 
corr_mtx(:,size(pvects,3)) = cm;

%computes spatial distance between all pixels
dists = nan(n_px^2, 1); 
count = 0;
for pix = 1:length(dists)
    count = count+1;
    [pix1_x, pix1_y] = ind2sub([sides sides], pixel_ids(pix,1));
    [pix2_x, pix2_y] = ind2sub([sides sides], pixel_ids(pix,2));
    dists(count) = pdist([pix1_x pix1_y; pix2_x, pix2_y]);
end

%sort similarity by dists
[dists, sort_idx] = sort(dists);
corr_mtx = corr_mtx(sort_idx, :);


%plot the average ctx
figure; hold on
avg_corm = mean(corr_mtx,2);
%plot(dists, avg_corm, 'o', 'color', [.8 .8 .8])

%draw line through means
u_dists = unique(dists(dists<pdist([0 0; sides/2 sides/2])));
dist_means = nan(length(u_dists), 1);
dist_stds = nan(length(u_dists), 1);
for id = 1:length(u_dists)
    dist_means(id) = mean(avg_corm(dists==u_dists(id)));
    dist_stds(id) = std(avg_corm(dists==u_dists(id)));
end
plot(u_dists, smooth(dist_means, 5), 'color', [0 .5 .75], 'linewidth', 3)
plot(u_dists, smooth(dist_means+dist_stds, 5), 'color', [0 .5 .75])
plot(u_dists, smooth(dist_means-dist_stds, 5), 'color', [0 .5 .75])

%}
end

%calculate similarity
function [corr_mtx_local, pixels_local] = simvect(pvects, ctxt, n_px)
    
    %trim and zscore pvects
    pvects = pvects(:,:,ctxt);
    pvects = zscore_mtx(pvects);

    %preallocate
    corr_mtx_local = nan(n_px^2,1); %r-vals
    pixels_local = nan(n_px^2,2); %col, row
    
    %pairwise correlations
    count_local=0;
    for i1 = 1:n_px
        for i2 = 1:n_px
            count_local = count_local+1;
            corr_mtx_local(count_local) = corr(pvects(i1, :)', pvects(i2, :)');
            pixels_local(count_local, :) = [i1, i2];
        end
    end
end
