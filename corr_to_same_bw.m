function corrs_out = corr_to_same_bw(bins)
%calculate population distance to same context type (e.g., b to b) as a function of visit
%number

%bins
%bins = 20;

%expected contexts
contexts_bw = [1 2 3 4];% 9 10 11 12 13];

%minimum confidence
min_conf = 3;

%regions
regions = [1 2];
%regions = 2;

%preallocate
all_vects_1a_w = [];
all_vects_1b_w = [];
all_vects_2a_w = [];
all_vects_2b_w = [];
all_vects_3a_w = [];
all_vects_3b_w = [];
all_vects_1a_b = [];
all_vects_1b_b = [];
all_vects_2a_b = [];
all_vects_2b_b = [];
all_vects_3a_b = [];
all_vects_3b_b = [];
corrs_out = nan(bins^2, 6);

%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%iterate through subjects
length_subjects = size(file_names_subjects{:},1);
for subject = 1:length_subjects
    
    %iterate through task folders
    rat = file_names_subjects{:}(subject,1).name
    task_folders = {'black_and_white' 'object_arrangement'};
    for task = 1 %bw only
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));

        %iterate through sessions
        length_sessions = size(file_list_sessions,1);
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
            
            %load session
            current_file = strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file));
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)), 'clusters', 'context_ids', 'eptrials')
            
            %skip sessions without cells
            if isempty(clusters)
                continue
            else
                %constrain clusters
                cluster_idx = clusters(:,2)>=min_conf & ismember(clusters(:,4), regions);
                clusters = clusters(cluster_idx, :);
                %skip sessions without cells
                if isempty(clusters)
                    continue
                end   
            end  
            
            %check if all task contexts were visited
            visited_contexts = unique(eptrials(~isnan(eptrials(:,6)),6));
            if ~isequal(visited_contexts(1:4), contexts_bw')
                continue
            end

            %get spike counts
            %contexts output as 4 cells in context_ids order
            context_ids = context_ids';
            %[spike_count] = get_allcounts(eptrials, clusters, context_ids, window_duration, min_conf);
            [~, ~, vectors_all] = pop_pixlecorr(eptrials, context_ids', clusters(:,1), bins);

            %load based on position distance between alike contexts 
            %(WITHIN CONTEXTS)
            %
            switch 1
                case ismember(context_ids, [1 2 3 4; 3 4 1 2], 'rows') %AABB
  
                    all_vects_1a_w = [all_vects_1a_w vectors_all(:,:,1)]; %first A
                    all_vects_1b_w = [all_vects_1b_w vectors_all(:,:,2)]; %second A
                    
                    all_vects_1a_w = [all_vects_1a_w vectors_all(:,:,3)]; %first B
                    all_vects_1b_w = [all_vects_1b_w vectors_all(:,:,4)]; %second B
                
                case ismember(context_ids, [1 3 2 4; 3 1 4 2], 'rows') %ABAB
                    
                    all_vects_2a_w = [all_vects_2a_w vectors_all(:,:,1)]; %first A
                    all_vects_2b_w = [all_vects_2b_w vectors_all(:,:,3)]; %second A
                    
                    all_vects_2a_w = [all_vects_2a_w vectors_all(:,:,2)]; %first B
                    all_vects_2b_w = [all_vects_2b_w vectors_all(:,:,4)]; %second B
                
                case ismember(context_ids, [1 3 4 2; 3 1 2 4], 'rows') %ABBA
                    
                    all_vects_3a_w = [all_vects_3a_w vectors_all(:,:,1)]; %first A
                    all_vects_3b_w = [all_vects_3b_w vectors_all(:,:,4)]; %second A
                    
                    all_vects_1a_w = [all_vects_1a_w vectors_all(:,:,2)]; %first B
                    all_vects_1b_w = [all_vects_1b_w vectors_all(:,:,3)]; %second B
                
            end
            %}
            
            %load based on position distance between distinct contexts
            %(BETWEEN CONTEXTS)
            %
            switch 1
                case ismember(context_ids, [1 2 3 4; 3 4 1 2], 'rows') %AABB
                
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,2)]; %second A
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,3)]; %first B
                    
                    all_vects_2a_b = [all_vects_2a_b vectors_all(:,:,1)]; %first A
                    all_vects_2b_b = [all_vects_2b_b vectors_all(:,:,3)]; %first B
                    
                    all_vects_2a_b = [all_vects_2a_b vectors_all(:,:,2)]; %second A
                    all_vects_2b_b = [all_vects_2b_b vectors_all(:,:,4)]; %second B
                    
                    all_vects_3a_b = [all_vects_3a_b vectors_all(:,:,1)]; %first A
                    all_vects_3b_b = [all_vects_3b_b vectors_all(:,:,4)]; %second B
 
                case ismember(context_ids, [1 3 2 4; 3 1 4 2], 'rows') %ABAB
                    
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,1)]; %first A
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,2)]; %first B
                    
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,2)]; %first B
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,3)]; %second A
                    
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,3)]; %second A
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,4)]; %second B
                    
                    all_vects_3a_b = [all_vects_3a_b vectors_all(:,:,1)]; %first A
                    all_vects_3b_b = [all_vects_3b_b vectors_all(:,:,4)]; %second B
                
                case ismember(context_ids, [1 3 4 2; 3 1 2 4], 'rows') %ABBA
                    
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,1)]; %first A
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,2)]; %first B
                    
                    all_vects_1a_b = [all_vects_1a_b vectors_all(:,:,3)]; %second B
                    all_vects_1b_b = [all_vects_1b_b vectors_all(:,:,4)]; %second A
                    
                    all_vects_2a_b = [all_vects_2a_b vectors_all(:,:,1)]; %first A
                    all_vects_2b_b = [all_vects_2b_b vectors_all(:,:,3)]; %second B
                    
                    all_vects_2a_b = [all_vects_2a_b vectors_all(:,:,2)]; %first B
                    all_vects_2b_b = [all_vects_2b_b vectors_all(:,:,4)]; %second A
                    
            end
            %}
            
        end
    end 
end

%COMPUTE CORRELATIONS
%

%standardize
all_vects_1a_w = zscore_mtx(all_vects_1a_w);
all_vects_1b_w = zscore_mtx(all_vects_1b_w);
all_vects_2a_w = zscore_mtx(all_vects_2a_w);
all_vects_2b_w = zscore_mtx(all_vects_2b_w);
all_vects_3a_w = zscore_mtx(all_vects_3a_w);
all_vects_3b_w = zscore_mtx(all_vects_3b_w);

all_vects_1a_b = zscore_mtx(all_vects_1a_b);
all_vects_1b_b = zscore_mtx(all_vects_1b_b);
all_vects_2a_b = zscore_mtx(all_vects_2a_b);
all_vects_2b_b = zscore_mtx(all_vects_2b_b);
all_vects_3a_b = zscore_mtx(all_vects_3a_b);
all_vects_3b_b = zscore_mtx(all_vects_3b_b);

%correlate
corrs_out(:,1) = popcorr(all_vects_1a_w, all_vects_1b_w);
corrs_out(:,3) = popcorr(all_vects_2a_w, all_vects_2b_w);
corrs_out(:,5) = popcorr(all_vects_3a_w, all_vects_3b_w);

%correlate
corrs_out(:,2) = popcorr(all_vects_1a_b, all_vects_1b_b);
corrs_out(:,4) = popcorr(all_vects_2a_b, all_vects_2b_b);
corrs_out(:,6) = popcorr(all_vects_3a_b, all_vects_3b_b);


%PLOT
%
corr_out_plot = cell(1, size(corrs_out,2));
for ip = 1:size(corrs_out,2)
    corr_out_plot{ip} = corrs_out(:,ip);
end

errorbar_barplot( corr_out_plot )
xticks(1:size(corrs_out,2))
xticklabels({'1w', '1b', '2w', '2b', '3w', '3b'})
xlabel('Lag Distance')
ylabel('Spatial Correlation')


%two way anova
%{
anova_groups1 = [ones(size(corrs_out(:,1))); repmat(2, size(corrs_out(:,2)));...
    ones(size(corrs_out(:,3))); repmat(2, size(corrs_out(:,4))); ...
    ones(size(corrs_out(:,5))); repmat(2, size(corrs_out(:,6)))]; %within/between

anova_groups2 = [ones(size(corrs_out(:,1))); ones(size(corrs_out(:,2)));...
    repmat(2, size(corrs_out(:,3))); repmat(2, size(corrs_out(:,4))); ...
    repmat(3, size(corrs_out(:,5))); repmat(3, size(corrs_out(:,6)))]; %lag position

[p, tbl] = anovan(corrs_out(:), {anova_groups1, anova_groups2}, 'varnames', {'within/between', 'lag position'})
%}


%INTERNAL FUNCTIONS
%

%pop distance function
%
    function c_out = popcorr(m1, m2)
        c_out = nan(size(m1,1),1);
        for si = 1:size(m1,1)
            c_out(si) = corr(m1(si,:)', m2(si,:)'); 
        end 
    end

end