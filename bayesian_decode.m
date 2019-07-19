function [class, posterior, f_x] = bayesian_decode(sample, training, group, tao, varargin)
% Bayesian decoder for spatial coding based on Kechen Zhang, Iris
% Ginzburg, Bruce L. McNaughton, Terrence J. Sejnowski, 1998.
%
% [CLASS, POSTERIOR] = BAYESIAN_DECODE(SAMPLE,TRAINING,GROUP, TAO, BINS) 
% classifies each row of the data in SAMPLE into one of the groups in TRAINING.  
% SAMPLE and TRAINING must be matrices with the same number of columns.  
% GROUP is a grouping variable for TRAINING.  Its unique values define 
% groups, and each element defines which group the corresponding row of 
% TRAINING belongs to. TAO is the size of the time window in seconds.
%
% Optional input BINS is the number of bins used, and will instigate 
%   Gaussian smoothing of the tuning curves, fi(x).
%
%   In practice:
%       'sample' and 'training' are matrices with time samples for rows and
%           cell firing rates for columns.
%       'train_ID' is a vector of pixel ids describing the rat was at each 
%           corresponding time sample (row) of the matrix 'training'.
%
% CLASS
%  indicates which group each row of SAMPLE has been assigned to, and is
%  of the same type as GROUP.
%
% POSTERIOR, a matrix
%  containing estimates of the posterior probabilities that the j'th
%  training group was the source of the i'th sample observation, i.e.
%  Pr{group j | obs i}.
%
%   In practice:
%       class is a vector of absolute guesses (in the form of pixel ids) as
%           to where the rat was at each time sample (row) of the matrix 
%           'sample'.
%       posterior is a matrix of distributed guesses (in the form of a
%           probability distribution over all pixel ids) to where the rat 
%           was at each time sample (row) of the matrix 'sample'. Rows
%           correspond to time samples, columns correspond to pixel ids.
%
%
%
%                                                       A.M.P. Miller, 2014
%


%classify inputs
%
    if nargin > 4
        switch nargin
            case 5
                bins = varargin{1};
            otherwise
                error(message('too many inputs')) 
        end
    end

    
%check inputs
%
    if size(group,1) ~= size(training,1)
        error('stats:classify:InputSizeMismatch',...
            'The length of GROUP must equal the number of rows in TRAINING.');
    elseif size(sample,2) ~= size(training,2)
        error('stats:classify:InputSizeMismatch',...
            'SAMPLE and TRAINING must have the same number of columns.');
    end

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% FROM THE TRAINING DATA %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Here, we use the training data to calculate the prior probability of the 
%   rat occupying each pixel, and the firing rate maps of each cell over 
%   all pixels.
%


    %PRIOR PROBABILITY, P(A), i.e. the spatial occupancy distribution 
    %
    % P(A) = N(A) / sum(N(A))
    %
    % N(A) is the number of timesamples with the rat in pixel A
    % 
    % This is the number of timesamples with the rat in each pixel A 
    % divided by the total number of timesamples
    %  
    % This gives us a row vector length pixels
    %

        %visted pixels
        [gindex,~,glevels] = grp2idx(group);
        ngroups = length(glevels);
        ncells = size(training,2); %used later
        ntrain = size(training,1); %used later
        nsamp = size(sample,1); %used later

        %how many times each pixel was visited
        gsize = hist(gindex,1:ngroups);

        %spatial occupancy
        P_A = (gsize / sum(gsize));
        
        
        %uniform prior
        P_A = ones(size(P_A));

    
    %FIRING RATE MAPS, f(x), i.e. tuning curve distribution
    %
    % f(x) is the mean firing rate observed in each pixel A for each cell
    %
    % This gives us a matrix pixles*cells
    %

        %preallocate
        f_x = nan(ngroups, ncells);

        %training
        for pixel = 1:ngroups
            f_x(pixel, :) = mean(training(gindex==pixel, :))/tao;
        end

        % P(B) can be smoothed with a Gaussian kernel
        if exist('bins', 'var')
            
            %Gaussian kernal mask, approx. 4cm^2
            krnl=gausswin(3)*gausswin(3)';
            mask=krnl./sum(sum(krnl));

            %prepare to reshape tuning curve into spatial position
            grid = reshape((1:(bins^2))', bins, bins);
            grid(~ismember(grid, glevels')) = NaN;
                
            %for size fix below
            size_fix = grid; 
            size_fix(~isnan(size_fix))=1;

            %smooth each cell independently
            for cell = 1:size(training,2)

                %reshape tuning curve into spatial position
                grid(glevels) = f_x(:,cell);
                spatial_cell_tuning = grid;
                
                %smooth
                %{
                smooth_iterations = 5; %CHANGE to repeat smooth
                for i = 1:smooth_iterations
                    spatial_cell_tuning = conv2nan(spatial_cell_tuning, mask, 'same');
                end
                spatial_cell_tuning = spatial_cell_tuning.*size_fix;%delete bloat
                %}
                %
                %spatial_cell_tuning = smooth2a(spatial_cell_tuning, 3);
                %}
                %reshape tuning curve back to original shape
                f_x(:, cell) = spatial_cell_tuning(~isnan(spatial_cell_tuning));
                
            end
        end

        %Poisson distributions (below) do not work with an expected value
        %of zero. Therefore, we will replace mean firing rates of 0 with an
        %arbitrary non-zero value, 1*10^-10 or 0.000000001.
        f_x(f_x==0) = 0.000000001;

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% FOR THE SAMPLE DATA %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Here, we use the above calculated distributions to calculate the
%   conditional, P(B|A), and then the posterior, P(A|B), for the sample 
%   data.
%
    
    %CONDITIONAL PROBABILITY, P(B|A), i.e. the likelihood
    %
    % This is the probability for the firing rates B to occur given that 
    %   the animal is in pixel A. If we assume that the rates have 
    %  	Poisson distributions (are as random as possible in time) and 
    %  	that different cells are statistically independent of one another, 
    %  	then we can obtain the explicit expression:
    %
    %           P(B|A) = multiply_over P(Bi|A)
    %               and
    %           P(Bi|A) = multiply_over (((tao*fi(A))^Bi)/Bi!)*e^(-tao*fi(A))
    %
    % Here, P(B) is not explicitly calculated, but is obtained by
    % 	assuming that P(A|B) must equal 1.
    %
    % This gives us a matrix sample_timesamples*pixels
    %

        %preallocate
        P_B_A_cell = nan(nsamp, ngroups, ncells); %each time window

        %expected number of spikes (tao * average firing of each
        %cell in each pixel)
        e_spikes = f_x*tao;
        
        %for each timestamp, calculate the conditional probability at each pixel for each cell
        %
        %may be able to cut 10s by using a 3d matrix instead of a loop
        %
        for tsamp = 1:nsamp

            %observed number of spikes (at this time sample, resized to 
            %match e_spikes)
            o_spikes = repmat(sample(tsamp, :), ngroups, 1);
            
            %fill conditional matrix
            %
            %This computes whole matrices at once for efficiency.
            %Transposed to make multiplication (below) easier.
            %
            %P_B_A_cell_p1(tsamp, pixel, :) = ((e_spikes^o_spikes)/factorial(o_spikes))*exp(-e_spikes);
            P_B_A_cell(tsamp, :, :) = exp(o_spikes.*log(e_spikes)...
                -e_spikes-gammaln(o_spikes+ones(size(o_spikes)))); %equivalent and more stable
            
        end

        %multiply across cells
        P_B_A = prod(P_B_A_cell, 3);

    %posterior probability, P(A|B)
    %
    % The reconstruction is based on the standard formula of conditional 
    %	probability: P(A|B) = P(B|A)*P(A) / P(B)
    %
    % Merging the above distributions gives us the formula for the posterior:
    %
    %       P(A|B) = P_A * P_B_A_cell * C(tao,B)
    %
    % Where C(tao,B) is a normilazation factor, such that sum(P(A|B)) is equal to 1.
    %
    % This gives us a matrix sample_timesamples*pixles
    
        %replicate rate map for each time sample
        P_A = repmat(P_A, nsamp, 1);

        %multiply at each timestamp and pixel
        P_A_B = P_A .* P_B_A;


        %apply correction, C(tao,B)
        for tsamp = 1:nsamp
            P_A_B(tsamp, :) = P_A_B(tsamp, :)./sum(P_A_B(tsamp, :));
        end

    %outputs
    %
    posterior = P_A_B;
    [~, class] = max(P_A_B, [], 2);%highest probability pixel
    origin_pixels = unique(group);
    class = origin_pixels(class);
    
end
    
    


