function [eptrials] = trials_pcia(event, pos, flags, context_ids, CSC)
%eptrials merges event and pos, and adds some additional columns:
%
%   eptrials(:,1) is a common timestamp in the form of 1:.01:N seconds with
%       additional timestamps corresponding to events of different kinds
%   eptrials(:,2) is the rat's X coordinate position data
%   eptrials(:,3) is the rat's Y coordinate position data
%   eptrials(:,4) is the spike event. Unique numbers correspond to unique
%       cells. Ones correspond to (non-spike) videosamples
%   eptrials(:,5)  contains the context visit number
%       1 = first context
%       2 = second context
%       3 = third context
%       4 = fourth context
%   eptrials(:,6) contains current context.
%       1 = black, first visit
%       2 = black, second visit
%       3 = white, first visit
%       4 = white, second visit
%       5 = arrangment 1
%       6 = arrangment 2
%       7 = object 1
%       8 = object 2
%       9 = pre-session
%       10 = iti1
%       11 = iti2
%       12 = iti3
%       13 = post-session
%   eptrials(:,7) contains the rat's velocity.
%   eptrials(:,8) contains the rat's acceleration
%   eptrials(:,9) contains the rat's head direction.
%   eptrials(:,10) contains proximity to object 1 (top, left)
%   eptrials(:,11) contains proximity to object 2 
%   eptrials(:,12) contains proximity to object 3 
%   eptrials(:,13) contains proximity to object 4 (bottom, right)
%   eptrials(:,14) contains head direction relative to object 1 
%       furthest left, upper object is at 0 degrees
%   eptrials(:,15) contains head direction relative to object 2 
%       middle left, lower object is at 0 degrees
%   eptrials(:,16) contains head direction relative to object 3 
%       middle right, upper object is at 0 degrees
%   eptrials(:,17) contains head direction relative to object 4 
%       furthest right, lower object is at 0 degrees
%   eptrials(:,18) contains TT 1.0 LFP samples (downsampled in GetCSC). THIS
%       WILL ONLY BE LOADED IF USER SPECIFIES SO IN 'loadnl.m'
%   eptrials(:,19) contains TT 9.0 LFP samples (downsampled in GetCSC). THIS
%       WILL ONLY BE LOADED IF USER SPECIFIES SO IN 'loadnl.m'


%MERGE event, pos, flags, and possibly CSC, matrices.
%
    if exist('CSC','var') %if CSC is loaded
        %expand to accommodate csc input
        event(:,end+1) = nan;
        pos(:,end+1) = nan;
        flags(:,end+1) = nan;
        
        %merge
        eptrials = sortrows([event;pos;flags;CSC], 1);
        
        %adds empty columns
        eptrials = [eptrials(:,1:5) ones(length(eptrials(:,1)),1) NaN(length(eptrials(:,1)),2) eptrials(:,6) NaN(length(eptrials(:,1)),3) eptrials(:,7)];
 
    else %if CSC is not loaded
        %merge
        eptrials = sortrows([event;pos;flags], 1);

        %adds empty columns
        eptrials = [eptrials(:,1:5) ones(length(eptrials(:,1)),1) NaN(length(eptrials(:,1)),2) eptrials(:,6)];
    end

    
%DELETE IMPOSSIBLE POSITIONS
%
eptrials(eptrials(:,2)>800 | eptrials(:,2)<50 | eptrials(:,3)>600 | eptrials(:,3) <1, 2) = NaN;
eptrials(eptrials(:,2)>800 | eptrials(:,2)<50 | eptrials(:,3)>600 | eptrials(:,3) <1, 3) = NaN; 
eptrials = interp_at_nans(eptrials);
    function eptrials = interp_at_nans(eptrials)
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

%DELETE POSITIONLESS ROWS
%
    eptrials(isnan(eptrials(:,2)) | isnan(eptrials(:,3)), :) = [];

%SMOOTH POSITION DATA
%
    eptrials(:,2) = smooth(eptrials(:,2),10);
    eptrials(:,3) = smooth(eptrials(:,3),10);
   
%CONTEXT VISIT NUMBER (column 5)
%
eptrials = trl_num(eptrials);
    function eptrials = trl_num(eptrials)
    % trl_num fills column 5 with the context visit number. 
    %
    % INPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
        
        %catch missed flag
        %
        eptrials = catch_missed_flag(eptrials);
        function eptrials = catch_missed_flag(eptrials) 
        %return an error if there are less than 7 or more than 8 flags. If
        %there are 7, add a flag at the most likely time.
        
            if sum(eptrials(:,5))==7
                figure; plot(eptrials(:,1), eptrials(:,5)); hold on

                    %identify session flag pairs (begining and end)
                    flag_times = eptrials(eptrials(:,5)==1,1);

                    %time between each flag
                    interflag_durations = diff(flag_times);

                    %too-long session
                    false_session = find(interflag_durations>(14*60));

                    %if start flag was missed
                    if interflag_durations(false_session-1) > (11.5*60)

                        %add flag 12min before end flag
                        new_flag_time = flag_times(false_session+1)-(12*60);
                        find_nearby_time = abs(eptrials(:,1) - repmat(new_flag_time, size(eptrials(:,1))));
                        nearby_time_index = eptrials(find(find_nearby_time==min(find_nearby_time),1),1);
                        eptrials(eptrials(:,1)==nearby_time_index,5)=1;
                        clearvars find_nearby_time

                        %report
                        warning('event flag added. see figure.') 

                    %if end flag was missed
                    elseif interflag_durations(false_session+1) > (11.5*60)

                        %add flag 12min after start flag
                        new_flag_time = flag_times(false_session)+(12*60);
                        find_nearby_time = abs(eptrials(:,1) - repmat(new_flag_time, size(eptrials(:,1))));
                        nearby_time_index = eptrials(find(find_nearby_time==min(find_nearby_time),1),1);
                        eptrials(eptrials(:,1)==nearby_time_index,5)=1;
                        clearvars find_nearby_time

                        %report
                        warning('event flag added. see figure.')

                    else
                        error('too few event flags. could not resolve.')
                    end

                plot([new_flag_time new_flag_time], [0 1], 'r-')
                title('new flag in red'); 
                xlabel('seconds');set(gca,'Ytick',[]);set(gca,'TickLength',[0, 0]);
                hold off

            elseif sum(eptrials(:,5))<7
                number_of_flags = nansum(eptrials(:,5))
                figure; plot(eptrials(:,1), eptrials(:,5))
                error('incorrect number of flags')
            end
        end
    
        %renumber flag intervals (context visits)
        %
        eptrials = renumber_trials(eptrials);
        function eptrials = renumber_trials(eptrials)
            flag_counter = 1;
            for i1 = 1:length(eptrials(:,5))
                %if there is a flag
                if eptrials(i1,5)==1
                    %increase flag interval by 1
                    flag_counter = flag_counter + 1;
                end
                %set trial number in column 5
                eptrials(i1,5) = flag_counter;
            end

            %delete periods of improbable lengths (not context visits)
            %{
            for trial = unique(eptrials(:,5))'
                starttime = min(eptrials(eptrials(:,5)==trial, 1));
                enddtime = max(eptrials(eptrials(:,5)==trial, 1));
                duration = enddtime - starttime;

                if duration < 11.5*60 %duration in seconds
                    eptrials(eptrials(:,5)==trial, 6) = 0;
                elseif trial == 1
                    eptrials(eptrials(:,5)==trial, 6) = 0;
                end
            end 
            %}
            
            eptrials(ismember(eptrials(:,5), [1 3 5 7 9]), 6) = 0;
            
            %renumber remaining trials (iti = 0; contexts = 1:4)
            eptrials(eptrials(:,6)==0, 5) = 0;
            flag_counter1 = 1;
            for trial1 = unique(eptrials(eptrials(:,6)==1,5))'
                %renumber
                eptrials(eptrials(:,5)==trial1 & eptrials(:,6)==1, 5) = flag_counter1;
                %increase trial number by 1
                flag_counter1 = flag_counter1 + 1;
            end
        end
        
        unique(eptrials(:,5))'
        
        %catch if too many context visits
        %
        eptrials = catch_extra_contexts(eptrials);
        function eptrials = catch_extra_contexts(eptrials)
            ctxt_nums = unique(eptrials(~isnan(eptrials(:,5)),5));
            ctxt_nums = ctxt_nums(ctxt_nums>0);
            if length(ctxt_nums) > 4

                %we can do something if there is one too many
                if length(ctxt_nums) == 5

                    %find starts and ends
                    ctxt_bounds = nan(length(ctxt_nums), 2);
                    for ictx = 1:length(ctxt_nums)
                        ctxt_bounds(ictx,:) = [min(eptrials(eptrials(:,5)==ctxt_nums(ictx), 1)) max(eptrials(eptrials(:,5)==ctxt_nums(ictx), 1))];
                    end
                    ctxt_bounds = ctxt_bounds(:);
                    ctxt_bounds_rounded = floor(ctxt_bounds);

                    %find repeat times
                    repeats = nan(1,2);
                    count = 1;

                    for icb = 1:length(unique(ctxt_bounds_rounded))
                       if sum(ctxt_bounds_rounded==ctxt_bounds_rounded(icb))>1
                           repeats(count) = ctxt_bounds(icb);
                           ctxt_bounds_rounded(icb) = 0;
                           count = count + 1;
                       end
                    end
                    repeats = sort(repeats);

                    %set eptrials(:,5) within bounds to 0
                    eptrials(eptrials(:,1)>=repeats(1) & eptrials(:,1)<=repeats(2), 5) = 0;

                    %renumber column 5
                    flag_count = 1;
                    for trial2 = unique(eptrials(eptrials(:,5)>0,5))'
                        %renumber
                        eptrials(eptrials(:,5)==trial2, 5) = flag_count;
                        %increase trial number by 1
                        flag_count = flag_count + 1;
                    end

                    %renumber column 6
                    eptrials(eptrials(:,5)==0,6) = 0;
                    eptrials(eptrials(:,5)>0,6) = 1;

                else
                    figure; plot(eptrials(:,5))
                    error('too many context visits')
                end

                %check that we didnt screw things up
                if length(unique(eptrials(~isnan(eptrials(:,5)) & eptrials(:,5)>0, 5))) ~= 4
                    figure; plot(eptrials(:,5))
                    error('incorrect number of context visits')
                end
            end
        end
    end


%CONTEXT ID (column 6)
%
eptrials = ctxt_id(eptrials, context_ids);
    function eptrials = ctxt_id(eptrials, context_ids)
    % ctxt_id fills column 6 with the context IDs. 
    %
    % INPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %

    
        %iti order is always the same; context order (context_ids) is input
        iti_ids = 9:13;
    
        %use counters to index iti and context id vectors
        iti_counter = 1;
        ctx_counter = 0;
        
        %used to identify changes between iti and contexts
        counter_hold = 0;
        
        %for all rows
        for i = 1:length(eptrials(:,6))
            
            %if there is a flag
            if eptrials(i,6)==0
                
                %number itis in order
                if counter_hold == 1
                    iti_counter = iti_counter + 1;
                    counter_hold = 0;
                end

                eptrials(i,6) = iti_ids(iti_counter);
                
            elseif eptrials(i,6)==1
                
                %number contexts according to the input context order
                if counter_hold == 0
                    ctx_counter = ctx_counter + 1;
                    counter_hold = 1;
                end
                
                eptrials(i,6) = context_ids(ctx_counter);
            end
        end
    end

  
%delete positions outside of contexts
eptrials(eptrials(:,5)<1,2:3) = nan;

%DELETE TRACKING JUMPS
%
for ctxt_visit = 1:4
	eptrials(eptrials(:,5)==ctxt_visit, 1:3) = delete_jumps(eptrials(eptrials(:,5)==ctxt_visit, 1:3), 1000, 5);
end
    function eptrials = delete_jumps(eptrials, too_big_origin, cum_mod)
    % delete_jumps attempts to correct position data by removing spurious
    % points creates by light interference. It also resamples the data at
    % 100hz.
    %
    % INPUT VARIABLES
    %   eptrials = the nueralynx video data with the first three columns 
    %      %time, xpos, and ypos
    %   too_big_origin = the max acceptable position change between adjacent
    %      %points. Larger position changes are deleted as noise.
    %   cum_mod = a correction factor that prevents deletions from creating
    %      %unacceptable position changes between points
    %
    % OUTPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %
    
        %prep to remove improbable changes in position
        too_big = too_big_origin;
        guilty_by_association = 4;

        %adjacent pos points for evaluating velocity
        p1 = eptrials(1, 2:3);
        p2 = eptrials(2, 2:3);

        %adjacent time points for evaluating velocity
        t1 = eptrials(1,1);
        t2 = eptrials(2,1);

        %preallocate
        deletions = zeros(length(eptrials(:,2)),1);
        dists = zeros(length(eptrials(:,2)),1);

        %iterate through adjacent points to evaluate velocity
        count = 0;
        for i = 1:length(eptrials(:,2))-2
            
            %velocity
            current_distance = pdist([p1; p2])/(t2-t1);
            dists(i+1) = current_distance;

            %if the current velocity is too big
            if current_distance > too_big

                %note that point (and the next 4) should be deleted (index for later)
                if length(eptrials(:,2))-2-i > guilty_by_association
                    deletions(i:i+guilty_by_association) = 1;
                end

                %move to the next point, but keep the first of the adjacent pair
                p2 = eptrials(i+2, 2:3);
                t2 = eptrials(i+2, 1);

                %each time it's too big, increase what is considered "too big"
                count = count + cum_mod;
                too_big = too_big + count;

            %if it's not too big
            else

                %reset what is considered "too big"
                too_big = too_big_origin;
                count = 0;

                %update points
                p1 = eptrials(i+1, 2:3);
                p2 = eptrials(i+2, 2:3);
                t1 = eptrials(i+1, 1);
                t2 = eptrials(i+2, 1);

            end
        end

        %index to delete dubious points
        deletions = logical(deletions);
        eptrials(deletions, 2:3) = NaN;

        %replace deleted points with interpolated values
        non_nan_pos = ~isnan(eptrials(:,2)) & ~isnan(eptrials(:,3)); %index
        
        known_time = eptrials(non_nan_pos, 1);
        known_place_2 = eptrials(non_nan_pos, 2);
        known_place_3 = eptrials(non_nan_pos, 3);
        
        [known_time_unq, kt_indexA] = unique(known_time);

        known_place_2_unq = known_place_2(kt_indexA);
        known_place_3_unq = known_place_3(kt_indexA);

        eptrials(:,2) = interp1(known_time_unq, known_place_2_unq, eptrials(:,1), 'linear');
        eptrials(:,3) = interp1(known_time_unq, known_place_3_unq, eptrials(:,1), 'linear');
        
    end

%SQUARE AND REPOSITION LOCATIONS
%{
eptrials = sqr_pos(eptrials);
    function eptrials = sqr_pos(eptrials)
        %normalizes the height and width of the (1mX1m) arena
        
        %for each context visit
        for visit = 1:4

            %axis bins
            bins = 10;
            xaxis_bins = linspace(floor(min(eptrials(eptrials(:,5)==visit,2))),ceil(max(eptrials(eptrials(:,5)==visit,2))),bins);
            yaxis_bins = linspace(floor(min(eptrials(eptrials(:,5)==visit,3))),ceil(max(eptrials(eptrials(:,5)==visit,3))),bins);

            %mins and maxes
            minmax_xy = nan(4,bins);
            
            %for each bin (in both axes)
            for bin = 1:(bins-1)

                %axis indices
                x_index = eptrials(:,2)>xaxis_bins(bin) & eptrials(:,2)<xaxis_bins(bin+1);
                y_index = eptrials(:,3)>yaxis_bins(bin) & eptrials(:,3)<yaxis_bins(bin+1);
                                
                %find min and max for x and y
                minmax_xy(1,bin) = min(eptrials(eptrials(:,5)==visit & y_index, 2));
                minmax_xy(2,bin) = max(eptrials(eptrials(:,5)==visit & y_index, 2));
                minmax_xy(3,bin) = min(eptrials(eptrials(:,5)==visit & x_index, 3));
                minmax_xy(4,bin) = max(eptrials(eptrials(:,5)==visit & x_index, 3));
            end 
            
                %sort and keep the fourth most extreme value
                extrm_val = 4;
                minmax_xy = sort(minmax_xy,2);
                minmax_xy = [minmax_xy(1,extrm_val) minmax_xy(2,bins-extrm_val); minmax_xy(3,extrm_val) minmax_xy(4,bins-extrm_val)];

                    %HARDCODED ROTATION
                    %
                    %corrects for the consistent misplacement of the arena
                    
                    %counterclockwise rotation angle in degrees
                    cc_deg = .2;
                    
                    %translate to radians
                    rotang = cc_deg/57.2957795;

                    %build rotation matrix (some calc II thing)
                    romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];

                    %rotation occurs around origin, so we temporarily center all points around the origin
                    eptrials(eptrials(:,5)==visit,2) = eptrials(eptrials(:,5)==visit,2) - mean(minmax_xy(1,:));
                    eptrials(eptrials(:,5)==visit,3) = eptrials(eptrials(:,5)==visit,3) - mean(minmax_xy(2,:));

                    %apply rotation
                    eptrials(eptrials(:,5)==visit,[2 3]) = (romat*eptrials(eptrials(:,5)==visit,[2 3])')';

                    %undo centering
                    eptrials(eptrials(:,5)==visit,2) = eptrials(eptrials(:,5)==visit,2) + mean(minmax_xy(1,:));
                    eptrials(eptrials(:,5)==visit,3) = eptrials(eptrials(:,5)==visit,3) + mean(minmax_xy(2,:));
                    
                %subtract mins from all position points
                eptrials(eptrials(:,5)==visit,2) = eptrials(eptrials(:,5)==visit,2) - repmat(minmax_xy(1,1), size(eptrials(eptrials(:,5)==visit,2)));
                eptrials(eptrials(:,5)==visit,3) = eptrials(eptrials(:,5)==visit,3) - repmat(minmax_xy(2,1), size(eptrials(eptrials(:,5)==visit,3)));
                minmax_xy = minmax_xy - [minmax_xy(1,1) minmax_xy(1,1); minmax_xy(2,1) minmax_xy(2,1)];
                
                %divide all points by their maximum (normalizing to 1)
                eptrials(eptrials(:,5)==visit,2) = eptrials(eptrials(:,5)==visit,2)./minmax_xy(1,2);
                eptrials(eptrials(:,5)==visit,3) = eptrials(eptrials(:,5)==visit,3)./minmax_xy(2,2);
                
                %delete out-of-bounds points
                eptrials((eptrials(:,2)<0 | eptrials(:,2)>1) & eptrials(:,5)==visit, 2) = nan;
                eptrials((eptrials(:,3)<0 | eptrials(:,3)>1) & eptrials(:,5)==visit, 3) = nan;
                
                %interp deleted points
                %
                %index for missing position values
                real_x = eptrials(~isnan(eptrials(:,2)) & eptrials(:,5)==visit,2);
                real_y = eptrials(~isnan(eptrials(:,3)) & eptrials(:,5)==visit,3);
                all_time = eptrials(eptrials(:,5)==visit,1);

                %only including unique time points among the missing rows
                [real_time_x, idx_x] = unique(eptrials(~isnan(eptrials(:,2)) & eptrials(:,5)==visit,1));
                [real_time_y, idx_y] = unique(eptrials(~isnan(eptrials(:,3)) & eptrials(:,5)==visit,1));

                %index real_x and real_y to match real_time unique elements
                eptrials(eptrials(:,5)==visit,2) = interp1(real_time_x, real_x(idx_x), all_time);
                eptrials(eptrials(:,5)==visit,3) = interp1(real_time_y, real_y(idx_y), all_time); 
        end
    end
%}

%VELOCITY (column 7)
%
[eptrials, span_vel] = vid_velocity(eptrials);
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
            eptrials(eptrials(:,5)==visit,7) = interp1(eptrials(original_vel_idx,1)', eptrials(original_vel_idx,7)', eptrials(eptrials(:,5)==visit,1)');
        
        end
    end

%ACCELERATION (column 8)
%
eptrials = vid_accleration(eptrials, span_vel);
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
            eptrials(eptrials(:,5)==visit,8) = interp1(eptrials(original_accl_idx,1)', eptrials(original_accl_idx,8)', eptrials(eptrials(:,5)==visit,1)');            
        end
    end


%HEAD DIRECTION (column 9)
%
eptrials = vid_hd(eptrials);
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
        for ctxt = 1:4;
            %existing HDs for this context
            original_HD_idx = ~isnan(eptrials(:,9)) & eptrials(:,5)==ctxt;
            %smooth original HDs
            eptrials(original_HD_idx,9) = circ_smooth(eptrials(original_HD_idx,9), [0 360], 40);
            %interp to all other time points
            eptrials(eptrials(:,5)==ctxt,9) = circ_interp1(eptrials(original_HD_idx,9)', [0 360], eptrials(eptrials(:,5)==ctxt,1)', eptrials(original_HD_idx,1)');
        end
    end


%                                  %
% CODE SPECIFIC TO OBJECT CONTEXTS %
%                                  %
%                                  %

%if it was an object-context day...
if sum(ismember(eptrials(:,6), [5 6 7 8]))>0
    
    %hard-coded object radius
    object_radius = .17/2; %meters
    
    %hard-coded average object locations for all sessions
    object_locations(1,:,1) = [mean([0.5 0.6875]) mean([0.6875 1-0.125])] + [-0.005 0];
    object_locations(2,:,1) = [mean([0.3125 0.5]) mean([0.5 0.6875])] + [-0.005 0.005];
    object_locations(3,:,1) = [mean([0.6875 1-0.125]) mean([0.3125 0.5])] + [0 0.005];
    object_locations(4,:,1) = [mean([0.3125 0.5]) mean([0.125 0.3125])] + [-0.005 0.005];
    object_locations(1,:,2) = [mean([0.5 0.6875]) mean([0.5 0.6875])] + [-0.0025 0.0025];
    object_locations(2,:,2) = [mean([0.6875 1-0.125]) mean([0.5 0.6875])] + [0 0.005];
    object_locations(3,:,2) = [mean([0.125 0.3125]) mean([0.3125 0.5])] + [-0.0075 0.005];
    object_locations(4,:,2) = [mean([0.3125 0.5]) mean([0.3125 0.5])] + [0 0.005];
    object_locations(1,:,3) = [mean([0.125 0.3125]) mean([0.6875 1-0.125])] + [-0.0025 0.005];
    object_locations(2,:,3) = [mean([0.6875 1-0.125]) mean([0.6875 1-0.125])] + [-0.0075 0.0075];
    object_locations(3,:,3) = [mean([0.125 0.3125]) mean([0.125 0.3125])] + [-0.0075 0.0025];
    object_locations(4,:,3) = [mean([0.6875 1-0.125]) mean([0.125 0.3125])] + [-0.005 -0.005];
    object_locations(1,:,4) = [mean([0.125 0.3125]) mean([0.6875 1-0.125])] + [0 0.005];
    object_locations(2,:,4) = [mean([0.6875 1-0.125]) mean([0.6875 1-0.125])] + [-0.0075 0.0075];
    object_locations(3,:,4) = [mean([0.125 0.3125]) mean([0.125 0.3125])] + [-0.005 0];
    object_locations(4,:,4) = [mean([0.6875 1-0.125]) mean([0.125 0.3125])] + [-0.005 -0.005];
    
    %OBJECT PROXIMITY AND HD ANGLE (columns 10:13)
    %calculates the distance from each rat position to each object position
    %for each context, cols 10:13, and the egocentric angle to each object,
    %cols 14:17
    %    
        %external index for brevity
        real_pos_idx = ~isnan(eptrials(:,2)) & ~isnan(eptrials(:,3));
        
        %iterate through ALL contexts
        for context = unique(eptrials(:,6))'

            %external index for brevity
            context_idx = eptrials(:,6)==context;
            
            %only object contexts
            if ismember(context, 5:8)
                
                %all rat positions in this context
                rat_pos = eptrials(real_pos_idx & context_idx, 2:3);
                
                %all rat positions in this context
                object_pos = object_locations(:,:,context-4)';
                
                %fill eptrials(10:13) with distance to center of objects
                eptrials(real_pos_idx & context_idx, 10:13) = dist(rat_pos, object_pos);

                %make it distance to edge of object by subtracting radius
                %negative values indicate climbing on object
                eptrials(real_pos_idx & context_idx, 10:13) = eptrials(real_pos_idx & context_idx, 10:13) - repmat(object_radius, size(eptrials(real_pos_idx & context_idx, 10:13)));
            
                %fill eptrials(14:17) with HD difference to center of objects
                pos_difference = repmat(reshape(object_pos,1,2,4), size(rat_pos,1), 1, 1) - repmat(rat_pos, 1, 1, 4); %x and y offset between rat and object
                eptrials(real_pos_idx & context_idx, 14:17) = atan2d(squeeze(pos_difference(:, 1, :)), squeeze(pos_difference(:, 2, :))); %angle between rat and object
                eptrials(real_pos_idx & context_idx, 14:17) = eptrials(real_pos_idx & context_idx, 14:17) + 360.*double(eptrials(real_pos_idx & context_idx, 14:17)<0); %force 360 degree scale
                eptrials(real_pos_idx & context_idx, 14:17) = eptrials(real_pos_idx & context_idx, 14:17) - repmat(eptrials(real_pos_idx & context_idx, 9), 1, 4); %subtract current HD
                eptrials(real_pos_idx & context_idx, 14:17) = eptrials(real_pos_idx & context_idx, 14:17) + 360.*double(eptrials(real_pos_idx & context_idx, 14:17)<0); %[0 360] add 360 to neg angles
        
            %non-object contexts get nans
            else
                eptrials(context_idx, 10:17) = nan;
            end
        end
end



%HOUSE CLEANING
%

%sets first time bin to 0.0000
eptrials(:,1) = eptrials(:,1) - ones(length(eptrials(:,1)),1).*eptrials(1,1);

end



        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

