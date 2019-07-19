function ctl_col_process_neuralynx
% Create a neurodata_ctl folder, with >rat>session sub-folders, then fill 
% subfolders with day#.mat files containing variables: file, eptrials,
% clusters. The day#.mat files are generated by processing neurodata folder
% files.


%path to neurodata_ctl (destination)
path_destination = '/Users/ampm/Documents/MATLAB/PCia/';

%get all the things in neurodata folder...
file_list_subjects = dir('neurodata_vc/');

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
    rat = file_names_subjects{:}(subject,1).name;

    %task folders
    task_folders = {'black_and_white' 'object_arrangement'};

    %iterate through tasks
    for task = 1:2

        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('neurodata_vc/', num2str(file_names_subjects{:}(subject,1).name),'/', task_folders{task}, '/*.mat'));

        %number of folders
        length_sessions = size(file_list_sessions,1);

        %iterate through sessions
        for session = 1:length_sessions

            day_file = file_list_sessions(session,1).name;
        
            
                %check if destination has a neurodata_ctl folder, if not, create it.
                if ~exist(strcat(path_destination, 'neurodata_ctl'), 'dir')
                    mkdir(strcat(path_destination, 'neurodata_ctl'))
                end

                %check if destination has corresponding rat folder, if not, create it.
                if ~exist(strcat(path_destination, 'neurodata_ctl/', rat), 'dir')
                    mkdir(strcat(path_destination, 'neurodata_ctl/', rat))
                end

                %check if destination has corresponding stage folder, if not, create it.
                if ~exist(strcat(path_destination, 'neurodata_ctl/', rat, '/', task_folders{task}), 'dir')
                    mkdir(strcat(path_destination, 'neurodata_ctl/', rat, '/', task_folders{task}))
                end

                %
                %check if destination has corresponding data file, if so, continue (skip).
                if exist(strcat(path_destination, 'neurodata_ctl/', rat, '/', task_folders{task}, '/', num2str(day_file)), 'file')
                    continue
                end
                %}

                %load session
                strcat(path_destination, 'neurodata_vc/', rat, '/', task_folders{task}, '/', num2str(day_file))
                load(strcat(path_destination, 'neurodata_vc/', rat, '/', task_folders{task}, '/', num2str(day_file)))

                    %screen files without cells
                    if ~isempty(clusters)
                        %remove affects of col var from eptrials
                        %   eptrials(:,7) contains the rat's velocity.
                        %   eptrials(:,8) contains the rat's acceleration
                        eptrials = ctrl_FR_by_col(eptrials, clusters(:,1), 8, 30, 500, 500);

                        %interp
                        for context_visit = 1:4
                            %pos
                            eptrials(eptrials(:,5)==context_visit, :) = interp_pos_nans(eptrials(eptrials(:,5)==context_visit, :));
                        end
                            %velocity
                            [eptrials, span_vel] = vid_velocity(eptrials);
                            %acceleration
                            eptrials = vid_accleration(eptrials, span_vel);
                            clear span_vel
                            %head direction
                            eptrials = vid_hd(eptrials);
                    end

                %save session in correct destination folder
                save(strcat(path_destination, 'neurodata_ctl/', rat, '/', task_folders{task}, '/', num2str(day_file)), 'eptrials', 'clusters', 'context_ids', 'origin_file')

                %clear variables
                clear eptrials
                clear clusters
                clear stem_runs
                clear clusts
                clear origin_file
                clear data
                clear context_ids
                
                
        end
        
    end

end


%interp functions
function eptrials = interp_pos_nans(eptrials)
%interpolate_at_nans finds missing position information in eptrials
%and interpolates to fill
%

    %index for missing position values
    real_x = eptrials(~isnan(eptrials(:,2)),2);
    real_y = eptrials(~isnan(eptrials(:,3)),3);
    all_time = eptrials(:,1);

    %only including unique time points among the missing rows
    [real_time_x, idx_x] = unique(eptrials(~isnan(eptrials(:,2)),1));
    [real_time_y, idx_y] = unique(eptrials(~isnan(eptrials(:,3)),1));

    %index real_x and real_y to match real_time unique elements
    eptrials(:,2) = interp1(real_time_x, real_x(idx_x), all_time);
    eptrials(:,3) = interp1(real_time_y, real_y(idx_y), all_time);  

end

function [eptrials, span_vel] = vid_velocity(eptrials)
%vid_velocity adds a column with the estimated instantaneous
%velocity at every video sample time-point (row). Other rows are nan.
%
% INPUTS
%   eptrials = the overall output; a matrix of all relevant information
%   cm_per_matlab_unit = conversion between xy coordniate space and real-life distance
%
% OUTPUT
%   eptrials = the overall output; a matrix of all relevant information
%    

%how many .01s vid frames to use for each average velocity
span_vel = 5;

%after normalization, 1 unit is 1m
m_per_matlab_unit = 1;

    %for each context
    for visit = 1:4

        %calculate distances between every video sample
        eptrials_vid = eptrials(eptrials(:,4)==1 & eptrials(:,5)==visit, 2:3);

        pos1 = [eptrials_vid(1:end-span_vel,1)'; eptrials_vid(1:end-span_vel,2)'];
        pos2 = [eptrials_vid((span_vel+1):end,1)'; eptrials_vid((span_vel+1):end,2)'];

        distances = nan(length(eptrials_vid),1);
        for i = 1:length(eptrials_vid)-span_vel
            distance = dist([pos1(:,i) pos2(:,i)]);
            distances(i+floor(span_vel/2)) = distance(1,2);
        end

        %calulate velocity from dists_overall (see above)
        velocity = (distances.*m_per_matlab_unit)./(0.01*span_vel);

        %replace nans at ends
        velocity(1:floor(span_vel/2)) = repmat(velocity(floor(span_vel/2)+1), size(velocity(1:floor(span_vel/2))));
        velocity((end-(floor(span_vel/2)+1)):end) = repmat(velocity(end-(floor(span_vel/2)+2)), size(velocity((end-(floor(span_vel/2)+1)):end)));

        %set rows
        eptrials(eptrials(:,4)==1 & eptrials(:,5)==visit,7) = smooth(velocity, 40);
        
        %interpolate nans
        original_vel_idx = ~isnan(eptrials(:,7)) & eptrials(:,5)==visit;
        real_v = eptrials(original_vel_idx,7);
        [real_time, idx_v] = unique(eptrials(original_vel_idx,1));
        eptrials(eptrials(:,5)==visit,7) = interp1(real_time, real_v(idx_v), eptrials(eptrials(:,5)==visit,1)');

    end
end

function eptrials = vid_accleration(eptrials, span_vel)
%vid_accleration adds a column with the estimated average
%acceleration at every video sample time-point (row). Other rows are nan.
%
% INPUTS
%   eptrials = the overall output; a matrix of all relevant information
%   span_vel = How many .01s vid frames were used to calculate each
%               average velocity in the vid_velocity function
%
% OUTPUT
%   eptrials = the overall output; a matrix of all relevant information
%    

%how many .01s vid frames to use for each avg accl
span_accl = 10;

    %for each context
    for visit = 1:4

        %calculate velocity differences between every video sample
        eptrials_vid = eptrials(eptrials(:,4)==1 & eptrials(:,5)==visit, 7);

        pos1 = eptrials_vid((floor(span_vel/2)+1):(end-span_accl))';
        %pos2 = eptrials_vid(((floor(span_vel/2)+1+1)+span_accl):end)';
        pos2 = eptrials_vid(span_accl:(end-(floor(span_vel/2)+1)))';

        acceleration = zeros(length(eptrials_vid),1);
        for i = 1:(length(eptrials_vid)-(floor(span_vel/2)+span_accl))
            difference = diff([pos1(:,i) pos2(:,i)]);
            acceleration(i+floor(span_accl/2)) = difference/(0.01*span_accl);
        end

        %set rows
        eptrials(eptrials(:,4)==1 & eptrials(:,5)==visit, 8) = smooth(acceleration, 10);

        %interpolate nans
        original_accl_idx = ~isnan(eptrials(:,8)) & eptrials(:,5)==visit;
        real_a = eptrials(original_accl_idx,8);
        [real_time, idx_a] = unique(eptrials(original_accl_idx,1));
        eptrials(eptrials(:,5)==visit,8) = interp1(real_time, real_a(idx_a), eptrials(eptrials(:,5)==visit,1)'); 
        
    end
end

function eptrials = vid_hd(eptrials)
%vid_hd adds a column with the estimated head direction.
%
% INPUTS
%   eptrials = the overall output; a matrix of all relevant information
%
% OUTPUT
%   eptrials = the overall output; a matrix of all relevant information
%        

    %only keep hd during context visits
    eptrials(eptrials(:,5)<1,9) = nan; 

    %smooth and interpolate within context visits
    %
    for ctxt = 1:4
        %existing HDs for this context
        original_HD_idx = ~isnan(eptrials(:,9)) & eptrials(:,5)==ctxt;
        %smooth original HDs
        eptrials(original_HD_idx,9) = circ_smooth(eptrials(original_HD_idx,9), [0 360], 40);
        %interp to all other time points
        eptrials(eptrials(:,5)==ctxt,9) = circ_interp1(eptrials(original_HD_idx,9)', [0 360], eptrials(eptrials(:,5)==ctxt,1)', eptrials(original_HD_idx,1)');
    end
end



end