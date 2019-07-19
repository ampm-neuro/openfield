function [H, contour_coords] = place_cell(eptrials, cluster, bins, figure_on, varargin)
%identifies place cells in a single session
%H is a vector of zeros and ones with the same length as clusts. H=1
%indicates that the corresponding cluster is a place field.

%INPUTS
%
%high or low input
%
    %0 for low, 1 for high
    if exist('varargin', 'var')
        hi_lo = varargin{1};
    else
        hi_lo = 1; %default to high
    end
    if ~ismember(hi_lo, [0 1])
        error('incorrect hi_lo input')
    end

    
%turns on process figures   
intermediate_plots = 0;

%RATE MATRIX
%
    %get firing rate matrix
    [~, spike_count, spatial_occupancy] = rate_mtx(eptrials, cluster, bins);

    %process matrix
    matrix = skagg_smooth(spike_count, spatial_occupancy); %smooth
    matrix = smooth2a(matrix,1); %smooth
    matrix = inpaint_nans(matrix); %interp nans

    %edges
    xedges = linspace(min(eptrials(:,2)), max(eptrials(:,2))+0.000000001, bins+1);
    yedges = linspace(min(eptrials(:,3)), max(eptrials(:,3))+0.000000001, bins+1);

    %process plot
    if intermediate_plots == 1
        figure; imagesc(matrix); axis square; axis off; title(['Cluster ' num2str(cluster)])
    end
    
    %visited pixles
    visited_pixles = length(matrix(matrix>0 & ~isnan(matrix)));


%DEFINING A PLACE FIELD
%
    %field firing rate
    rate_vector = sort(matrix(~isnan(matrix)));
    if hi_lo == 1
        field_rate = rate_vector(ceil(length(rate_vector)*.75));
    else
        field_rate = rate_vector(ceil(length(rate_vector)*.25));
    end

    %minimum size relative to visited area (contiguous pixles)
    min_in_sqr_cms = 36;
    pxl_sqr_cms = (100/bins)^2;
    min_size = min_in_sqr_cms / pxl_sqr_cms;

    %maximum size relative to visited area
    max_size = round(visited_pixles.*0.5);

    %within:without firing rate ratio (zscore diff between in and out
    %means)
    field_ratio = 2;

    %reliability (proportion of passes on which the cell fired above field_rate)
    reliability_min = 0.50;



%IDENTIFY PLACE FIELDS
%
    %RATE
    %
    %pass matrix through minimum rate conditions
    if hi_lo == 1
        logic_ratematrix = matrix > field_rate;
    else
        logic_ratematrix = matrix < field_rate;
    end
    %process figure
    if intermediate_plots == 1
        figure; imagesc(logic_ratematrix); title rate; axis square
    end

    
    %SIZE
    %
    %identify candidate fields (contiguous points)
    %calls the function 'contiguous'
    [logic_contiguitymatrix] = contiguous(double(logic_ratematrix), bins);
    
    %delete fields smaller than minimum size or larger than maximum size
    for field = 1:max(logic_contiguitymatrix(:))    
        if sum(sum(logic_contiguitymatrix == field)) < min_size
            logic_contiguitymatrix(logic_contiguitymatrix == field) = 0;
        elseif sum(sum((logic_contiguitymatrix == field))) > max_size
            logic_contiguitymatrix(logic_contiguitymatrix == field) = 0;
        end
    end
    remaining_fields = unique(logic_contiguitymatrix(logic_contiguitymatrix>0));
    
    %process figure
    if intermediate_plots == 1
        figure; imagesc(logic_contiguitymatrix); title size; axis square
    end
    
    
%WITHIN:WITHOUT RATE
%
    %delete fields with insufficient within:without firing rate
    for f = 1:length(remaining_fields)
        field = remaining_fields(f);

        infield_rate_mean = nanmean(matrix(logic_contiguitymatrix == field));
        %infield_rate_std = nanstd(matrix(logic_contiguitymatrix == field));
        outfield_rate_mean = nanmean(matrix(logic_contiguitymatrix ~= field));
        outfield_rate_std = nanstd(matrix(logic_contiguitymatrix ~= field));
        
        %how many standard deviations between in and out of field rates
        %zscore_diff = (infield_rate_mean-outfield_rate_mean) / mean([infield_rate_std outfield_rate_std])
        zscore_diff = (infield_rate_mean-outfield_rate_mean) / outfield_rate_std;

        if abs(zscore_diff) < field_ratio
            logic_contiguitymatrix(logic_contiguitymatrix == field) = 0;
        end
    end
    remaining_fields = unique(logic_contiguitymatrix(logic_contiguitymatrix>0));
    
    %process figure
    if intermediate_plots == 1
        figure; imagesc(logic_contiguitymatrix); title w/w-rate; axis square
    end
    
%RELIABILITY
%    
    %delete fields with insufficient reliability
    for f = 1:length(remaining_fields)
        field = remaining_fields(f);

        %time_binpos = [times binned_xpos binned_xpos event]
        [~, xbinz] = histc(eptrials(eptrials(:,4)==1 | eptrials(:,4)==cluster, 2), xedges);
        [~, ybinz] = histc(eptrials(eptrials(:,4)==1 | eptrials(:,4)==cluster, 3), yedges);
            xbinz(xbinz == 0) = 1;
            ybinz(ybinz == 0) = 1;
        time_binpos = [eptrials(eptrials(:,4)==1 | eptrials(:,4)==cluster, 1) xbinz ybinz eptrials(eptrials(:,4)==1 | eptrials(:,4)==cluster, 4)];

        %call function reliability
        [reliability_score] = reliability(logic_contiguitymatrix, time_binpos, field, cluster, field_rate, hi_lo);

        if reliability_score < reliability_min
            logic_contiguitymatrix(logic_contiguitymatrix == field) = 0;
        end
    end
    remaining_fields = unique(logic_contiguitymatrix(logic_contiguitymatrix>0));

    %process figure
    if intermediate_plots == 1
        figure; imagesc(logic_contiguitymatrix); title reliability; axis square
    end
    
    
%LABEL EACH FIELD
%
    %renumber fields AND COUNT AND PLOT
    for f = 1:length(remaining_fields)
        field = remaining_fields(f);
        %re-lable
        logic_contiguitymatrix(logic_contiguitymatrix == field) = f;
    end
    remaining_fields = unique(logic_contiguitymatrix(logic_contiguitymatrix>0));
    
    
%LOCATE FIELDS
%
    for field = 1:length(remaining_fields)
        %peaks
        rate_vect = sort(matrix(logic_contiguitymatrix == field));
        [row, col] = find(matrix == rate_vect(floor(.99*length(rate_vect))) & logic_contiguitymatrix == field);
    end
    
    %center of mass
    %XXX
    
    
%CALCULATE OUTPUT
%
    %output
    if ~isempty(remaining_fields)
        H = length(remaining_fields);
    else
        H = 0;
    end
    
    %preallocate additional outputs
    %{
    field_sizes = nan(length(remaining_fields), 1);
    rate_ratios = nan(length(remaining_fields), 1);
    
    %load additional outputs
    count = 0;
    for field = 1:length(remaining_fields)
        count = count+1;
        field_sizes(count) = sum(sum((logic_contiguitymatrix == field)));
        rate_ratios(count) = nanmean(matrix(logic_contiguitymatrix == field))/nanmean(matrix(logic_contiguitymatrix ~= field));
    end
    %}
 
    
%CONTOUR COORDINATES
%
    %find outline of fields 
    %(contour coordinates using contour function)
    contour_coords = cell(length(remaining_fields));
    for field = 1:length(remaining_fields)
        contour_coords{field} = contourc(padarray(double(logic_contiguitymatrix == field), [1 1], 'both'), 1); %pad bc contour hates edges
        contour_coords{field}(:,1) = [];
        end_field = find(contour_coords{field}(1,:)==.5, 1, 'first');
        contour_coords{field}(:, end_field:end) = [];
    end
    
%PLOT
%    
    %prep figure
    if ~isempty(remaining_fields) && figure_on > 0
        figure
        axis([0 1 0 1]); axis square
        set(gca,'TickLength',[0, 0]); box off
    end

    for field = 1:length(remaining_fields)

        %plot field countours (without heatmap)
        if figure_on == 1

            %convert heatmap coords to video coords
            for pxl = 1:size(contour_coords,2)
                if contour_coords(1,pxl) - round(contour_coords(1,pxl)) == 0
                    contour_coords(1,pxl) = xedges(contour_coords(1,pxl));
                else
                    contour_coords(1,pxl) = mean([xedges(floor(contour_coords(1,pxl))) xedges(ceil(contour_coords(1,pxl)))]);
                end

                if contour_coords(2,pxl) - round(contour_coords(2,pxl)) == 0
                    contour_coords(2,pxl) = yedges(contour_coords(2,pxl));
                else
                    contour_coords(2,pxl) = mean([yedges(floor(contour_coords(2,pxl))) yedges(ceil(contour_coords(2,pxl)))]);
                end
            end

            %remove padding from coords
            contour_coords = contour_coords - ones(size(contour_coords))./(bins);
            contour_coords = contour_coords + (ones(size(contour_coords))./(bins))/2;
            
            %smooth coords
            smooth_wndw = 8;
            contour_coords_sm(1,:) = smooth([contour_coords(1,:) contour_coords(1,1)],smooth_wndw);
            contour_coords_sm(2,:) = smooth([contour_coords(2,:) contour_coords(2,1)],smooth_wndw);

            %plot contour
            hold on; plot(contour_coords_sm(1,:), contour_coords_sm(2,:))

            %plot rate peak
            hold on; plot(xedges(col), yedges(row), '.', 'Color', [0 0 0], 'markersize', 25)
            
        end

        %plot field countours over heatmap
        if figure_on == 2

            if field == 1
               imagesc(matrix)
            end

            %remove padding from coords
            contour_coords = contour_coords - ones(size(contour_coords));
            
            %correct for bins
            contour_coords(1,:) = contour_coords(1,:) + repmat(.5, size(contour_coords(1,:)));
            contour_coords(2,:) = contour_coords(2,:) + repmat(.5, size(contour_coords(2,:)));
            
            %smooth coords
            smooth_wndw = 8;
            contour_coords_sm(1,:) = smooth([contour_coords(1,:) contour_coords(1,1)],smooth_wndw);
            contour_coords_sm(2,:) = smooth([contour_coords(2,:) contour_coords(2,1)],smooth_wndw);
            
            %plot contour
            hold on; plot(contour_coords_sm(1,:), contour_coords_sm(2,:), 'k-', 'linewidth', 2.5)

            %plot rate peak
            hold on; plot(col-.5, row+.5, '.', 'Color', [0 0 0], 'markersize', 25)

            title(['Place Field(s) for Cluster ' num2str(cluster)])
            axis([-1 bins+1 -1 bins+1])
            axis square
            axis off
            colorbar

            clear countour_coords_sm
            
        end
    end


    
%INTERNAL FUNCTIONS
%
    %rate matrix
        function [rate_matrix, spike_count, spatial_occupancy] = rate_mtx(eptrials, cluster, bins)
        %function [rate_matrix, spike_count, spatial_occupancy] = rate_mtx(eptrials, cluster, bins)
        %
        % rate_mtx is a simple heatmap plot with adaptive smoothing. It computes
        % the firing rate of the neuron 'cluster' within a grid with 'bins' rows and
        % 'bins' columns. It depends on precise formatting of the input eptrials,
        % where column 1 is time, column 2 is x position, column 3 is y position,
        % and column four contains cluster events labeled by their cluster ID, as
        % well as timestamps labeled as 1.
        %
        % rate_mtx outputs the raw rate matrix 'rate_matrix', which is computed by
        % divided the output 'spike_count' by the output 'spatial_occupancy'.
        %
        % rate_mtx also plots a heatmap of the smoothed rate matrix. It uses the
        % adaptive smoothing technique first employed by Skagg et al.
        %
        % ampm 2017


            %check inputs
            if ~ismember(cluster, unique(eptrials(:,4)))
                error('Cluster does not appear in eptrials.')
            end

            %all spike events in two vectors
            xs = eptrials(eptrials(:,4)==cluster, 2);
            ys = eptrials(eptrials(:,4)==cluster, 3);

            %all time samples in two vectors 
            xt = eptrials(eptrials(:,4)==1, 2);
            yt = eptrials(eptrials(:,4)==1, 3);

            %evenly spaced bins of x and y coordinate ranges
            pos_edges_x = linspace(min(eptrials(:,2)), max(eptrials(:,2))+0.000000001, bins+1);
            pos_edges_y = linspace(min(eptrials(:,3)), max(eptrials(:,3))+0.000000001, bins+1);

            %2d histogram of event counts
            spike_count = histcounts2(ys, xs, pos_edges_y, pos_edges_x);
            spatial_occupancy = histcounts2(yt, xt, pos_edges_y, pos_edges_x)./100;

            %flip
            spike_count = flipud(spike_count);
            spatial_occupancy = flipud(spatial_occupancy);

            %divide spikes by time for rate
            rate_matrix = spike_count./spatial_occupancy;
        end
    
    %smoothing functions
        function matrixOut = smooth2a(matrixIn,Nr,Nc)
        % Smooths 2D array data.  Ignores NaN's.
        %
        %function matrixOut = smooth2a(matrixIn,Nr,Nc)
        % 
        % This function smooths the data in matrixIn using a mean filter over a
        % rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
        % element "i" by the mean of the rectange centered on "i".  Any NaN
        % elements are ignored in the averaging.  If element "i" is a NaN, then it
        % will be preserved as NaN in the output.  At the edges of the matrix,
        % where you cannot build a full rectangle, as much of the rectangle that
        % fits on your matrix is used (similar to the default on Matlab's builtin
        % function "smooth").
        % 
        % "matrixIn": original matrix
        % "Nr": number of points used to smooth rows
        % "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
        % 
        % "matrixOut": smoothed version of original matrix
        % 
        % 
        % 	Written by Greg Reeves, March 2009.
        % 	Division of Biology
        % 	Caltech
        % 
        % 	Inspired by "smooth2", written by Kelly Hilands, October 2004
        % 	Applied Research Laboratory
        % 	Penn State University
        % 
        % 	Developed from code written by Olof Liungman, 1997
        % 	Dept. of Oceanography, Earth Sciences Centre
        % 	G?teborg University, Sweden
        % 	E-mail: olof.liungman@oce.gu.se

            %
            % Initial error statements and definitions
            %
            if nargin < 2, error('Not enough input arguments!'), end

            N(1) = Nr; 
            if nargin < 3, N(2) = N(1); else N(2) = Nc; end

            if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
            if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

            %
            % Building matrices that will compute running sums.  The left-matrix, eL,
            % smooths along the rows.  The right-matrix, eR, smooths along the
            % columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
            % (2*Nc+1) rectangle centered on element "i".
            %
            [row_,column] = size(matrixIn);
            eL = spdiags(ones(row_,2*N(1)+1),(-N(1):N(1)),row_,row_);
            eR = spdiags(ones(column,2*N(2)+1),(-N(2):N(2)),column,column);

            %
            % Setting all "NaN" elements of "matrixIn" to zero so that these will not
            % affect the summation.  (If this isn't done, any sum that includes a NaN
            % will also become NaN.)
            %
            A = isnan(matrixIn);
            matrixIn(A) = 0;

            %
            % For each element, we have to count how many non-NaN elements went into
            % the sums.  This is so we can divide by that number to get a mean.  We use
            % the same matrices to do this (ie, "eL" and "eR").
            %
            nrmlize = eL*(~A)*eR;
            nrmlize(A) = NaN;

            %
            % Actually taking the mean.
            %
            matrixOut = eL*matrixIn*eR;
            matrixOut = matrixOut./nrmlize;
        end
        function mtx_smth = skagg_smooth(spike_counts, occupancy_counts)
        % An adaptive smoothing method to optimize the trade-off between blurring 
        % error and sampling error. The firing rate at each bin was estimated
        % by expanding a circle around the point until the radius of the circle 
        % (in bins) is greater than a constant (try 10000) divided by n*sqrt(s),
        % where n is the number of occupancy samples within the circle and s is 
        % the total number of spikes within the circle. With a position sampling 
        % rate of 100 Hz, the firing rate at that point was then set to 100*n*s.
        %
        %Skaggs, McNaughton, Wilson, Barnes 1996, see also Ito 2015. This method is
        %favored by the Mosers.

        %skagg constant
        skg_c = 10000;

        occupancy_counts = occupancy_counts.*100;
        occupancy_counts(occupancy_counts==0) = nan;

        %preallocate rate matrix
        mtx_smth = nan(size(spike_counts));

            %iterate through rows
            for ir = 1:size(spike_counts,1)

                %iterate through columns
                for ic = 1:size(spike_counts,2)

                    %skip (leave nan) if no occupancy
                    if isnan(occupancy_counts(ir, ic))
                        continue
                    end

                    %preset test values counter
                    radius = 0;
                    skagg_val = 1; %arbitrarily higher than radius

                    %keep trying until pass skagg test
                    while radius < skagg_val

                        %expand smooth mask size
                        radius = radius+1;

                        %Circle with radius centered at ir,ic
                        [mgx, mgy] = meshgrid(1:size(spike_counts,1), 1:size(spike_counts,2));
                        circle_idx = sqrt((mgx-ic).^2+(mgy-ir).^2)<=radius;

                        %sum within circle area
                        circ_spikes = nansum(spike_counts(circle_idx));
                        circ_occupancy = nansum(occupancy_counts(circle_idx));

                        %calculate skagg test
                        skagg_val = skg_c / (circ_occupancy * sqrt(circ_spikes));

                    end

                    %set smoothed rate
                    mtx_smth(ir, ic) = 100 * circ_spikes / circ_occupancy;
                end

            end
        end

    %contiguous
        function [logic_contiguitymatrix] = contiguous(matrix, bins)
        %finds contiguous patches of ones that are atleast a certain size (defined
        %below)

            %SET MINIMUM SIZE HERE
            %field always requires a square of atleast ~4cm^2 pixles
            size_x = ceil((4*bins)/100);
            size_y = ceil((4*bins)/100);

            %start positions
            x=1;y=1;

            %FIRST PASS: find patches of 1s of minimum size
            %
            while 1

                %if every item in minimum size region is nonzero
                if isempty(find(matrix(y:y+size_y, x:x+size_x) == 0, 1))

                    %change items to 2
                    matrix(y:y+size_y, x:x+size_x) = 2;

                    x = x+1;
                else
                    x = x+1;
                end

                %at end of row, shift down one and begin again on left
                if x+size_x > length(matrix(1,:))
                    x = 1; y = y+1;
                end

                %at end of last row, break
                if y+size_y > length(matrix(:,1))
                    break
                end

            end

            %counter to see when secton pass is complete
            complete = 1;

            x=1;y=1;
            %SECOND PASS: include contiguous 1s
            %
            while 1

                %spotlight
                spotlight = matrix(y:y+1, x:x+1);

                %if 2x2 spotlight includes at least two 2s and at least one 1
                if sum(sum(spotlight == 2))>1 && ~isempty(find(spotlight == 1, 1))

                %change 1s to 2s
                    spotlight(spotlight == 1) = 2;
                    matrix(y:y+1, x:x+1) = spotlight;

                    x = x+1;

                    %reset completeness measure
                    complete = 1;

                else
                    x = x+1;

                    %increase completeness measure
                    complete = complete+1;

                    if complete > length(matrix(:))
                        break
                    end

                end

                %at end of row, shift down one and begin again on left
                if x+1 > length(matrix(1,:))
                    x = 1; y = y+1;
                end

                %at end of last row, break
                if y+1 > length(matrix(:,1))
                    x = 1; y = 1;
                end

            end

            %set to just 1's (field) and 0's (not field)
            matrix(matrix == 1) = 0;
            matrix(matrix == 2) = 1;

            %THIRD PASS: Give a unique number to each isolated field
            %
            %start positions
            id = 1;
            x=1;y=1;
            while 1

                %iterating a spotlight to identify pixles in need of changing
                if isequal(matrix(y:y+size_y, x:x+size_x), ones(size(matrix(y:y+size_y, x:x+size_x))))

                    %next unique field ID
                    id = id+1;

                    %Find first unIDed field coordinate (first 1)
                    yc = y;
                    xc = x;

                    %Change first 1 and all contiguous 1's to a unique, higher number
                    %flood_fill(xc,yc,xc,yc,id,1,size(fill,2),size(fill,1));
                    field_items = flood_fill(matrix, yc, xc);
                    matrix(field_items) = id;

                    x = x+1;
                else
                    x = x+1;
                end

                %at end of row, shift down one and begin again on left
                if x+size_x > length(matrix(1,:))
                    x = 1; y = y+1;
                end

                %at end of last row, break
                if y+size_y > length(matrix(:,1))
                    %matrix
                    break
                end

            end


            %subtrack all nonzero numbers by 1
            matrix = matrix - (matrix>0);


            %clear fill
            logic_contiguitymatrix = matrix;

        end

    %flood fill
        function ms=flood_fill(matrix,yc,xc)
        % flood fill, scan line algoritm
        % ms=flood_fill_new(I,r,tol)
        % ms - pixels numbers that flooded
        % numbering: number=y+matrix_size_Y*(x-1), x y - pixels coordinaties

            % (yc, xc)- first point of selection

            [matrix_size_Y, matrix_size_X] = size(matrix); % image size

            % stak, where seeds will be stored
            maximum_stack_size=10000; %needs to be higher for crazier shapes
            stak = zeros(maximum_stack_size, 2, 'int32');
            stak(1,1) = xc;
            stak(1,2) = yc; %first seed
            remaining_stack = 1; %starting stack length

            % establish pixel store
            maximum_field_size = 1000000; % margin
            field_pixles=zeros(maximum_field_size,1,'int32'); % predifined array to increase speed
            field_size=0; % 0 points initially

            % to pixel number %??? %number=y+matrix_size_Y*(x-1), x y - pixels coordinaties
            tn=@(xx,yy) yy+matrix_size_Y*(xx-1);


            while true
                % get seed from stack:
                x_seed = stak(remaining_stack,1); %x coord of seed
                y_seed = stak(remaining_stack,2); %y coord of seed
                remaining_stack = remaining_stack-1;


                % LINE SCAN TO RIGHT

                %These variables alternate below true and false below, and their
                %current state is used to gate code behavior
                sku=false; % seed key, true if seed added, up
                skd=false; % same for down
                sku1=false; % this key is need to prewnt extra seed when move to left, up
                skd1=false; % same for left

                %xtt - active pixle within the row
                for xtt = x_seed:matrix_size_X %from x coord of seed to end of line
                    if matrix(y_seed,xtt) == 1 %if current pixle needs to be changed
                        % add pixels to field
                        field_size=field_size+1;
                        field_pixles(field_size)=tn(xtt,y_seed);
                    else
                        break;
                    end

                    % TRY TO ADD (just one) SEED UP
                    if y_seed~=matrix_size_Y %if we're not already at the top
                        if matrix(y_seed+1,xtt) == 1
                            if ~sku %if we don't already have a seed from this line

                                if all(tn(xtt,y_seed+1)~=field_pixles(1:field_size)) % if free space (none of the above pixles are already in the field?)
                                    % add to stack
                                    remaining_stack=remaining_stack+1;
                                    stak(remaining_stack,1)=xtt;
                                    stak(remaining_stack,2)=y_seed+1;
                                    sku=true; %now we have a seed
                                end
                            end
                        else
                            sku=false; %don't have a seed
                        end
                        if xtt==x_seed
                            sku1=sku; % memorize (that we're not already at the top?), will be used when to left
                        end
                    end

                    % TRY TO ADD (just one) SEED DOWN
                    if y_seed~=1 %if we're not already at the bottom
                        if matrix(y_seed-1,xtt) == 1
                            if ~skd
                                if all(tn(xtt,y_seed-1)~=field_pixles(1:field_size)) % if free space
                                    % add to stack
                                    remaining_stack=remaining_stack+1;
                                    stak(remaining_stack,1)=xtt;
                                    stak(remaining_stack,2)=y_seed-1;
                                    skd=true;
                                end
                            end
                        else
                            skd=false;
                        end
                        if xtt==x_seed
                            skd1=skd; % memorize, will be used when to left
                        end
                    end
                end



                % LINE SCAN TO LEFT
                %sku=false; % seed key, true if seed added
                %skd=false;
                sku=sku1;
                skd=skd1;
                if x_seed~=1
                    for xtt=(x_seed-1):-1:1 

                        if matrix(y_seed,xtt) == 1 %if current pixle needs to be changed
                            % add pixel
                            field_size=field_size+1;
                            field_pixles(field_size)=tn(xtt,y_seed);
                        else
                            break;
                        end

                        % TRY TO ADD (just one) SEED UP
                        if y_seed~=matrix_size_Y
                            if matrix(y_seed+1,xtt) == 1
                                if ~sku
                                    if all(tn(xtt,y_seed+1)~=field_pixles(1:field_size)) % if free space
                                        % add to stack
                                        remaining_stack=remaining_stack+1;
                                        stak(remaining_stack,1)=xtt;
                                        stak(remaining_stack,2)=y_seed+1;
                                        sku=true;
                                    end
                                end
                            else
                                sku=false;
                            end
                        end

                        % TRY TO ADD (just one) SEED DOWN
                        if y_seed~=1
                            if matrix(y_seed-1,xtt) == 1
                                if ~skd
                                    if all(tn(xtt,y_seed-1)~=field_pixles(1:field_size)) % if free space
                                        % add to stack
                                        remaining_stack=remaining_stack+1;
                                        stak(remaining_stack,1)=xtt;
                                        stak(remaining_stack,2)=y_seed-1;
                                        skd=true;
                                    end
                                end
                            else
                                skd=false;
                            end
                        end
                    end
                end

                if remaining_stack==0 % no more seed
                    break; % stop
                end


            end

            %output
            ms=field_pixles(1:field_size);
        end

    %reliability
        function [reliability_score] = reliability(logic_contiguitymatrix, time_binpos, field, cluster, field_rate, hi_lo)
        %detirmines the proportion of passes through field during which the cell fires

            %counters
            success = 0;
            pass = 0;

            %time bounds
            pass_start = nan;
            pass_end = nan;

            %identify set of coordinates coorresponding to field
            [i,j] = find(logic_contiguitymatrix == field);

            %iterate through every time point
            for time = 1:length(time_binpos(:,1))

                %if the rat's position at the current time point is within the field
                if sum(time_binpos(time,2) == j & time_binpos(time,3) == i) > 0

                    %if we haven't found a pass start yet 
                    if isnan(pass_start)

                        %record start time
                        pass_start = time;

                        %erase previous end time
                        pass_end = nan;

                    end

                %if the rat is NOT in the field
                else

                    %if we haven't found pass end yet (but do have a start)
                    if ~isnan(pass_start) && isnan(pass_end)

                        %record end time
                        pass_end = time;

                        %increase pass counter
                        pass = pass + 1;

                        %firing rate during window pass_start : pass_end
                        starttime = time_binpos(pass_start,1);
                        endtime = time_binpos(pass_end,1);
                        time_idx = time_binpos(:,1) >= time_binpos(pass_start,1) & time_binpos(:,1) <= time_binpos(pass_end,1);
                        num_of_spikes = sum(time_binpos(:,4) == cluster & time_idx);
                        time_in_pass = endtime - starttime;
                        w_fr = num_of_spikes / time_in_pass;

                        if hi_lo == 1 && w_fr > field_rate
                            success = success +1;
                        elseif hi_lo == 0 && w_fr < field_rate
                            success = success +1;
                        end
                        

                        %reset pass_start
                        pass_start = nan;
                    end
                end
            end

            %reliability (high numbers = more reliability)
            %must have at least 5 passes!
            if pass < 5
                reliability_score = 0;
            else
                reliability_score = success/pass;
            end

        end
    

end