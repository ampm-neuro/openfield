function [p_correct, assignments, distances, class_test] = md_classify_traintest(data_train, class_train, data_test, class_test)
%performs a minimum distance classification of each row of data_mtx.
%reports success rate (p_correct) and the classified assingments
%(assigngments). class is list of correct classifications.

%dimensionality correction
dc = sqrt(size(data_train,2));

%remove classes that are not in training set
trained_classes = unique(class_train);

data_test = data_test(ismember(class_test, trained_classes), :);
class_test = class_test(ismember(class_test, trained_classes), :);

%preallocate
assignments = nan(size(class_test));
distances = nan(size(data_test,1), length(unique(class_test)));

%calculate class means from training set (training means)
training_means = nan(length(trained_classes), size(data_train,2));
for itm = trained_classes'
    training_means(trained_classes==itm, :) = mean(data_train(class_train==itm, :));
end

%iterate through rows of the test data
for ir = 1:size(data_test,1)
    
    %current sample
    cs = data_test(ir,:);
    
    %distance from sample to training means
    dists_hold = pdist([cs; training_means]);
    dists_hold = dists_hold(1:size(training_means,1));
    dists_hold = dists_hold./dc; %dim correction
    distances(ir, :) = dists_hold;
    
    %classify
    [~, min_idx] = min(dists_hold);
    assignments(ir) = trained_classes(min_idx);
    
end

%proportion correctly classified
p_correct = sum(assignments == class_test) / length(class_test);




end