function [rotations, new_corrs, old_corrs] = all_rot_fit(rate_dists)
%runs rot_fit on every neuron in the nested cell 'rate_dists'

%preallocate
rotations = nan(length(rate_dists),1);
new_corrs = nan(size(rotations));
old_corrs = nan(size(rotations));

    for coi=1:length(rate_dists) 

    %average black contexts and average white contexts
    %to make function input
    rd = [mean([rate_dists{coi}{1,1};rate_dists{coi}{1,2}]); ...
        mean([rate_dists{coi}{1,3};rate_dists{coi}{1,4}])];

    %run rotation function
    [rotation, rot_cor, orig_cor] = rot_fit(rd);

    if ~isempty(rotation)
        %load output
        rotations(coi) = rotation;
        new_corrs(coi) = rot_cor;
        old_corrs(coi) = orig_cor;
    end

end