function [p_correct, assignments, distances] = md_classify_inputmeans(data_mtx, class, tms)
%performs a minimum distance classification of each row of data_mtx.
%reports success rate (p_correct) and the classified assingments
%(assigngments). class is list of correct classifications.
%
%tms are training means

%dimensionality correction
dc = sqrt(size(data_mtx,2));

%preallocate
assignments = nan(size(class));
distances = nan(size(data_mtx,1), size(tms,1));
classes = 1:size(tms,1);


%iterate through rows of data_mtx (samples)
for ir = 1:size(data_mtx,1)
    
    %current sample
    cs = data_mtx(ir,:);
    
    %training set (remove current sample)
    data_mtx_c = data_mtx(setdiff(1:size(data_mtx,1), ir), :);
    class_c = class(setdiff(1:size(data_mtx,1), ir));
    
    %training means
    %input
    
    %distance from sample to training means
    dists_hold = pdist([cs; tms]); 
    dists_hold = dists_hold(1:size(tms,1));
    dists_hold = dists_hold./dc;
    distances(ir, :) = dists_hold;
    
    %classify
    [~, min_idx] = min(dists_hold);
    assignments(ir) = classes(min_idx);
    
end

%proportion correctly classified
p_correct = sum(assignments == class) / length(class);

%figure; imagesc([cs; tms])
%[cs; tms]

end