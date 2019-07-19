function [high_shuffle, low_shuffle, shuffle_distribution] = ...
    miller_infoscore_shuffle(ds, number_of_shuffles, pval)
%miller_infoscore(ds)
%
%see miller_infoscore. shuffles across all bins (within and between).
%pval is two sided
%

shuffle_distribution = nan(number_of_shuffles,1);
ds_vect = ds(:);
for ishuf = 1:length(shuffle_distribution)

    %shuffle
    ds_vect = ds_vect(randperm(length(ds_vect)));
    ds = reshape(ds_vect, size(ds));
    
    
    %var within
    var_within = nan(size(ds,1),1);
    for irow = 1:size(ds, 1)
        %excludes nans
        var_within(irow) = nanvar(ds(irow,:));
    end
    var_within = mean(var_within);


    %var between
    var_between = nan(size(ds,2),1);
    for icol = 1:size(ds, 2)
        %excludes nans
        var_between(icol) = nanvar(ds(:,icol));
    end
    var_between = mean(var_between);


    %calculate ratio
    shuffle_distribution(ishuf) = var_within/(var_within + var_between);

end
shuffle_distribution = sort(shuffle_distribution);

high_shuffle_idx = ceil((1-pval/2)*number_of_shuffles);
low_shuffle_idx = floor((pval/2)*number_of_shuffles);
if low_shuffle_idx == 0
    low_shuffle_idx = 1;
end
high_shuffle = shuffle_distribution(high_shuffle_idx);
low_shuffle = shuffle_distribution(low_shuffle_idx);


end