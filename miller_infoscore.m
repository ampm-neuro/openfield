function [mis, max_pairwise_r]  = miller_infoscore(ds)
%miller_infoscore(ds)
%
%ds is the distribution matrix, with rows corresponding to unique samples
%of the same distribution, and columns corresponding to unique points 
%within the distribution
%
%e.g., ds could contain firing rates at multiple spatial locations. in that
%case, each row would be from one recording session, and each column would
%correspond to a single spatial bin.
%
%measure of reliability and specificity.
%
%takes the ratio of the variance between points in a distribution and the
%variance between points across samples (var within)/(var between)
%
%high values indicate that there are unique values at different points
%within a sample, and that these same value are seen across multiple
%samples = high info score.
%

%standardize rates within rows
%ds = ds - repmat(min(ds, [], 2), 1, size(ds,2));
%ds = ds ./ repmat(max(ds, [], 2), 1, size(ds,2));

%var within
var_within = nan(size(ds,1),1);
for irow = 1:size(ds, 1)
    %excludes nans
    var_within(irow) = nanvar(ds(irow,:));
end
var_within = nanmean(var_within);


%var between
var_between = nan(size(ds,2),1);
for icol = 1:size(ds, 2)
    %excludes nans
    var_between(icol) = nanvar(ds(:,icol));
end
var_between = nanmean(var_between);


%calculate ratio
mis = var_within/(var_within + var_between);


%correlations between samples
pw_comps = unique(sort(pairwise(1:size(ds,1)),2), 'rows');
pairwise_rs = nan(size(pw_comps,1), 2);
for ipw = 1:size(pw_comps,1)
    nonan_idx = ~isnan(ds(pw_comps(ipw,1),:)) & ~isnan(ds(pw_comps(ipw,2),:));
    noinf_idx = ~isinf(ds(pw_comps(ipw,1),:)) & ~isinf(ds(pw_comps(ipw,2),:));
    [r,p] = corr(ds(pw_comps(ipw,1), nonan_idx & noinf_idx)', ds(pw_comps(ipw,2), nonan_idx & noinf_idx)');
    pairwise_rs(ipw, :) = [r p];
end
max_pairwise_r = [mean(pairwise_rs(:,1)) max(pairwise_rs(:,2))];


end