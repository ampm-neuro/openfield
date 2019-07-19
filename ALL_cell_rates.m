function [all_rates_bw, all_rates_first_bw, all_rates_second_bw, all_rates_ia, all_rates_first_ia, all_rates_second_ia] = ALL_cell_rates
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.


%, all_rates_itis, all_rates_first_itis, all_rates_second_itis

all_rates_bw = [];
all_rates_first_bw = [];
all_rates_second_bw = [];

all_rates_ia = [];
all_rates_first_ia = [];
all_rates_second_ia = [];


%{
all_rates_itis = [];
all_rates_first_itis = [];
all_rates_second_itis = [];
%}
%names%get all the things in neurodata folder...
file_list_subjects = dir('neurodata/');

%hard coded erasure of irrelevant directory folders
file_list_subjects(1:2) = [];

%exclude non-folders
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%number of folders
length_subjects = size(file_names_subjects{:},1);
%file_names_subjects{:}(1:length_subjects,1).name;

%iterate through subjects
for subject = 1:length_subjects

    %print update
    rat = file_names_subjects{:}(subject,1).name
    
    %get all the *.mat in subject folder...
    file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/*.mat'));

    
    %number of folders
    length_sessions = size(file_list_sessions,1);
    
        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name
            %day = session
            
            %load session
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/' ,num2str(day_file)));

            
            %clusters
            %clusters = clusters(clusters(:,2)==1, 1);
            clusters = clusters(:, 1);
            
            if isempty(clusters)
                continue
            end
            %{
            rate_check = zeros(size(clusters,1),1);
            for rc = 1:size(clusters,1)
               if  sum(eptrials(:,4)==clusters(rc),1)/max(eptrials(:,1)) < 3
                    rate_check(rc) = 1;
               end
            end
            
            clusters = clusters(logical(rate_check), :);
            
            
            if isempty(clusters)
                continue
            end
            %}
            
            %{
                %INSERT FUNCTION HERE
                [contexts, itis] = context_flags(eptrials);
                %itis_ordered = nan(size(itis));
                %for i = 1:length(ctxts_ordered)-1 
                %    itis_ordered(i) = min(itis(itis>ctxts_ordered(i)));
                %end
                %clearvars i;
                bwbw;

                sel_contexts = contexts(abab_order);

                [cell_rates, cell_rates_first, cell_rates_second] = ctxt_rates(eptrials, clusters, sel_contexts);
                %[cell_rates_itis, cell_rates_first_itis, cell_rates_second_itis] = ctxt_rates(eptrials, clusters, itis);

                if bwbw == 1
                    
                    %heatmaps_trial(eptrials, clusters(:,1),sel_contexts) 
                    
                    all_rates_bw = [all_rates_bw; cell_rates];
                    all_rates_first_bw = [all_rates_first_bw; cell_rates_first];
                    all_rates_second_bw = [all_rates_second_bw; cell_rates_second];
                elseif bwbw == 0
                    all_rates_ia = [all_rates_ia; cell_rates];
                    all_rates_first_ia = [all_rates_first_ia; cell_rates_first];
                    all_rates_second_ia = [all_rates_second_ia; cell_rates_second];
                end
                   
                %all_rates_first = [all_rates_first; cell_rates_first];
                %all_rates_second = [all_rates_second; cell_rates_second];

                %all_rates_itis = [all_rates_itis; cell_rates_itis];
                %all_rates_first_itis = [all_rates_first_itis; cell_rates_first_itis];
                %all_rates_second_itis = [all_rates_second_itis; cell_rates_second_itis];
            
            %}
            
            [num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell)]
            
        end
    
end


end

