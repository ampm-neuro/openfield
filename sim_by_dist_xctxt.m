function sim_by_dist_xctxt(pvects)
%load pop-pxlecorr.mat
%computes the similarity (rval) between pixles as a function of distance

%number of pixels in pvects
n_px = size(pvects,1);
sides = sqrt(n_px);

%subsample pvects
%{
lo_hi = [0.2 0.9];
pvects_idx = reshape(1:n_px, sides, sides);
pvects_idx = pvects_idx(round(lo_hi(1)*sides : lo_hi(2)*sides), round(lo_hi(1)*sides : lo_hi(2)*sides));
pvects = pvects(pvects_idx(:), :,:);
n_px = size(pvects,1);
sides = sqrt(n_px);
%}

%compute similarity between all pixels
%
%within context
withins = [1 2; 3 4];
corr_mtx_withins = nan(n_px^2, size(withins,1));
for comparison = 1:size(withins,1)
    ctxts = withins(comparison, :);
    [cm, pixel_ids] = simvect(pvects, ctxts, n_px); 
    corr_mtx_withins(:,comparison) = cm; 
end
%
%between context
betweens = [1 3; 1 4; 2 3; 2 4];
corr_mtx_betweens = nan(n_px^2, size(withins,1));
for comparison = 1:size(betweens,1)
    ctxts = betweens(comparison, :);
    corr_mtx_betweens(:,comparison) = simvect(pvects, ctxts, n_px); 
end


%computes spatial distance between all pixels
dists = nan(n_px^2, 1); 
count = 0;
for pix = 1:length(dists)
    count = count+1;
    [pix1_x, pix1_y] = ind2sub([sides sides], pixel_ids(pix,1));
    [pix2_x, pix2_y] = ind2sub([sides sides], pixel_ids(pix,2));
    dists(count) = pdist([pix1_x pix1_y; pix2_x, pix2_y]);
end


%sort similarity by dists and average
[dists, sort_idx] = sort(dists);
corr_mtx_withins = corr_mtx_withins(sort_idx, :);
corr_mtx_withins = mean(corr_mtx_withins,2);
corr_mtx_betweens = corr_mtx_betweens(sort_idx, :);
corr_mtx_betweens = mean(corr_mtx_betweens,2);

%plot the average ctx
figure; hold on
u_dists = unique(dists(dists<pdist([0 0; sides/2 sides/2])));
dist_means_withins = nan(length(u_dists), 1);
dist_stds_withins = nan(length(u_dists), 1);
dist_means_betweens = nan(length(u_dists), 1);
dist_stds_betweens = nan(length(u_dists), 1);
for id = 1:length(u_dists)
    dist_means_withins(id) = mean(corr_mtx_withins(dists==u_dists(id)));
    %dist_stds_withins(id) = std(corr_mtx_withins(dists==u_dists(id)));
    dist_stds_withins(id) = std(corr_mtx_withins(dists==u_dists(id)))./sqrt(sum(dists==u_dists(id)));
    dist_means_betweens(id) = mean(corr_mtx_betweens(dists==u_dists(id)));
    %dist_stds_betweens(id) = std(corr_mtx_betweens(dists==u_dists(id)));
    dist_stds_betweens(id) = std(corr_mtx_betweens(dists==u_dists(id)))./sqrt(sum(dists==u_dists(id)));
end

plot(u_dists, smooth(dist_means_withins, 3), 'color', [0 .5 .75], 'linewidth', 3)
plot(u_dists, smooth(dist_means_withins+dist_stds_withins, 3), 'color', [0 .5 .75])
plot(u_dists, smooth(dist_means_withins-dist_stds_withins, 3), 'color', [0 .5 .75])

plot(u_dists, smooth(dist_means_betweens, 3), 'color', [0.85 0.3 0.1], 'linewidth', 3)
plot(u_dists, smooth(dist_means_betweens+dist_stds_betweens, 3), 'color', [0.85 0.3 0.1])
plot(u_dists, smooth(dist_means_betweens-dist_stds_betweens, 3), 'color', [0.85 0.3 0.1])

%}
end

%calculate similarity
function [corr_mtx_local, pixels_local] = simvect(pvects, ctxt, n_px)
    
    %trim and zscore pvects
    pvects1 = pvects(:,:,ctxt(1));
    pvects1 = zscore_mtx(pvects1);
    pvects2 = pvects(:,:,ctxt(2));
    pvects2 = zscore_mtx(pvects2);

    %preallocate
    corr_mtx_local = nan(n_px^2,1); %r-vals
    pixels_local = nan(n_px^2,2); %col, row
    
    %pairwise correlations
    count_local=0;
    for i1 = 1:n_px
        for i2 = 1:n_px
            count_local = count_local+1;
            corr_mtx_local(count_local) = corr(pvects1(i1, :)', pvects2(i2, :)');
            pixels_local(count_local, :) = [i1, i2];
        end
    end
end
