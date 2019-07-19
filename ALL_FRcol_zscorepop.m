function [pvals, rvals, poly_vals] = ALL_FRcol_zscorepop
% Get all *.mat files out of each stage folder out of each subject folder
% out of neurodata.

count = 0;
pvals = [];
rvals = [];
poly_vals = [];

rates_all = [];
col_all = [];


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
    
    %task folders
    task_folders = {'black_and_white' 'object_arrangement'};
    
    for task = 1%:2
    
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));


        %number of folders
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions


            day_file = file_list_sessions(session,1).name
            %day = session
            
            %load session
            load(strcat('/Users/ampm/Documents/MATLAB/PCia/neurodata/',num2str(rat),'/', task_folders{task},'/' ,num2str(day_file)));
                if isempty(clusters)
                    continue
                end
            
            
            % CLUSTERS
            % clusters(:,1) contains cell id numbers tt.cluster
            % clusters(:,2) contains confidence 1-5
            %                   3 is an includable cell
            %                   4 is a very clear cluster
            %                   5 is a figure-worthy example cluster
            % clutsers(:,3) contains cell type based on waveform
            %                   1 is a pyramidal cell
            %                   2 is an interneuron
            %                   3 is unknown (probably inter)
            % clusters(:,4) contains RSC region during recording
            %                   1 is dysgranular (RSA)
            %                   2 is granular B (RGB)
            %
            clusters = clusters(clusters(:,2)>=3 & ismember(clusters(:,4), [1 2]), 1);
            if isempty(clusters)
                continue
            end
            
            
            
            %for all cells
            count = count+length(clusters);
            %rates_all = [];
            %col_all = [];
            ctxt_id = [];

            for ctxt = 1:4

                cur_eptrials = eptrials(eptrials(:,5)==ctxt, :);

                %   eptrials(:,7) contains the rat's velocity.
                %   eptrials(:,8) contains the rat's acceleration
                %   eptrials(:,9) contains the rat's head direction.
                col = 8;
                [rates_out, col_out] = FR_by_col(cur_eptrials, clusters, col, 30, 500, 500);
                
                %norm
                rates_out = zscore_mtx(rates_out);
                
                %remove mean
                %rates_out = rates_out - repmat(nanmean(rates_out), size(rates_out,1), 1); 
                
                
                
                %combine
                rates_all = [rates_all; rates_out(:)];
                col_all = [col_all; repmat(col_out, size(rates_out,2), 1)];
                ctxt_id = [ctxt_id; repmat(mode(cur_eptrials(:,6)), size(col_out))];

            end
                
            %for each cell  
            %{
            for c = 1:length(clusters)
                
                %correllation
                %[r, p, poly_val] = fit_line(col_all, rates_all(:,c), 1);
                [r, p, poly_val] = fit_line(col_all, rates_all(:,c), 2);
                pvals = [pvals; p];
                rvals = [rvals; r];     
                poly_vals = [poly_vals; poly_val];
                
                %plot by context BW
                %{
                hold on;
                for ictxt = unique(ctxt_id)'
                    ctxt_idx = ctxt_id == ictxt;
                   if ismember(ictxt, [1 2])
                       plot(col_all(ctxt_idx), rates_all(ctxt_idx,c), '.', 'markersize', 8, 'color', [0 0 0])
                   elseif ismember(ictxt, [3 4])
                       plot(col_all(ctxt_idx), rates_all(ctxt_idx,c), '.', 'markersize', 8, 'color', [.7 .7 .7])
                   end
                end
                %}
                
                box off
                set(gca,'TickLength',[0, 0]);
                title(['rat ' num2str(rat) ', session ' num2str(day_file) ', cell ' num2str(clusters(c))])
                %ylim([-3.5 3.5])
                ylabel('Firing Rate')
                if col == 7
                    xlabel('Velocity (m/s)')
                elseif col == 8
                    xlabel('Acceleration')
                elseif col == 9
                    xlabel('Head direaction')
                end
                    
                %title([num2str(rat) ' ' num2str(session) ' ' num2str(c) ' in context ' num2str(mode(cur_eptrials(:,6))) ' , visit ' num2str(ctxt) ])
                %saveas(gcf,['c' num2str(count)], 'epsc')
                %
                
            end
            %}
            
            
        end
    end
end

%correllation
%col_all = abs(col_all);
[r, p, poly_val] = fit_line(col_all, rates_all, 2);
pvals = [pvals; p];
rvals = [rvals; r];     
poly_vals = [poly_vals; poly_val];

%plot by context BW
%{
hold on;
for ictxt = unique(ctxt_id)'
    ctxt_idx = ctxt_id == ictxt;
   if ismember(ictxt, [1 2])
       plot(col_all(ctxt_idx), rates_all(ctxt_idx,c), '.', 'markersize', 8, 'color', [0 0 0])
   elseif ismember(ictxt, [3 4])
       plot(col_all(ctxt_idx), rates_all(ctxt_idx,c), '.', 'markersize', 8, 'color', [.7 .7 .7])
   end
end
%}

box off
set(gca,'TickLength',[0, 0]);
%title(['rat ' num2str(rat) ', session ' num2str(day_file) ', cell ' num2str(clusters(c))])
%ylim([-3.5 3.5])
ylabel('Firing Rate (zscored hz)')
if col == 7
    xlabel('Velocity (m/s)')
elseif col == 8
    xlabel('Acceleration')
elseif col == 9
    xlabel('Head direaction')
end


end

