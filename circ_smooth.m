function [smooth_vector, ideal_vector] = circ_smooth(varargin)
%
% function [smooth_vector, ideal_vector] = circ_smooth(rough_vector, min_max, smooth_window_size)
%
% circ_smooth applies the basic matlab smoothing function on circular data
%
% circ_smooth forgoes any knowledge about circular mathematics in favor of 
% brute force. Hopefully you have a shiny computer.
% 
% ALL INPUTS ARE ROW VECTORS
%
% circ_smooth assumes that the max value in 'rough_vector' is contiguous 
% with the minimum value (e.g., 360degrees is contiguous with 1degrees when
% calculating angles around a circle). If this is not true (e.g., there are
% no items in rough_vector with exactly 360 degrees, leaving the actual min_max 
% 1 degree short (0 - 359)), you can assert the min_max using the input value
% 'min_max' (e.g., [0 360]).
%
% The general strategy is to pad rough_vector above the max and below the min 
% with copies of itself. An ideal vector is then calculated by iterating through 
% the columns of the padded_matrix and selecting the item in the current column 
% with the minimum distance to the selected item in the previous column. Basic 
% smoothing function is then run on the ideal vector. Finally, we correct the 
% items in the ideal vector back to their original position.
%
%EXAMPLE DEMONSTRATING STRATEGY: 
%
%   rough_vector = [1 2 3 4 1 2 3]   
%   min_max = [1 4]
%
% Pad data
%   padded_matrix = [ 5  6  7  8  5  6  7; %in practice there are
%                     1  2  3  4  1  2  3;  usually more than 3 rows
%                    -3 -2 -1  0 -3 -2 -1]
%
% Calculate ideal_vector
%   ideal_vector = [1 2 3 4 5 6 7]
%
% Smooth ideal_vector
%   smooth_ideal_vector = smooth([1 2 3 4 5 6 7]) %applies basic smoothing
%
% Correcting back to original positions, add row to ideal_vector
%   ideal_vector = [1 2 3 4 5 6 7]
%
% Subtract appropriate 360degree shift
%   smooth_vector = [1 2 3 4 1 2 3] %correct back to original positions
%



%CHECK INPUTS
%

switch nargin

    case 0
        error(message('Need more input arguments'))
    case 1
        rough_vector = varargin{1};
        disp('min_max will be [min(rough_vector) max(rough_vector)]')
        min_max = [min(rough_vector) max(rough_vector)];
        smooth_window = 5;
    case 2
        rough_vector = varargin{1};
        min_max = varargin{2};
        
        
        %check for proper min_max input
        if ~isequal(size(min_max), [1 2])
            error(message('min_max must be a 1x2 row vector'));
        end
        smooth_window = 5;
        
    case 3
        rough_vector = varargin{1};
        min_max = varargin{2};
        
        %check for proper min_max input
        if ~isequal(size(min_max), [1 2])
            error(message('min_max must be a 1x2 row vector'));
        end
        
        smooth_window = varargin{3};
  
    otherwise
        error(message('Too many input arguments'))
end

%check for rough_vector orientation
if size(rough_vector,1) > size(rough_vector,2) 
    rough_vector = rough_vector';
end
if ~ismember(1, size(rough_vector))
    error(message('rough_vector must be a row vector'));
end

%SOME PRE-CODING THINGS
%

%pad_vector is a local function that calculates the pad vector by shifting 
%rough_vector up or down
function shifted_vector = shift_vector(vector, shift, min_max)

    min_max_shift = ones(1,length(vector)).*(min_max(2)*shift);
    shifted_vector = vector + min_max_shift;

end



%BUILD PADDED_MATRIX
%

%pads determines the number of rows that surround rough_vector to make 
%the padded_matrix. Must be an even number.
pads = 100;

%preallocate padded_matrix
padded_matrix = nan(pads+1, length(rough_vector));

%build padded_matrix with shifted vectors
for row = 1 : size(padded_matrix, 1)
    
    %pad_vector is a local function
    padded_matrix(row,:) = shift_vector(rough_vector, row - ceil(size(padded_matrix, 1)/2), min_max);

end



%CALCULATE IDEAL_VECTOR
%

%preallocate ideal_vector
ideal_vector = nan(1, length(rough_vector));

%start with the first item in rough_vector
ideal_vector(1,1) = rough_vector(1);

%find the next item
for current_item = 2:length(ideal_vector)
       
    %indexing out just the column under consideration
    current_column = padded_matrix(:, current_item);

  	%a column full of the current item
  	temp_item_column = ones(size(current_column)).*ideal_vector(1, current_item-1);
    
  	%how far is each item in the current column from the current item?
   	distances = abs(current_column - temp_item_column);

  	%fill ideal_vector with the closest item and the shift the corresponds 
   	%to the row it came from
        
    %skipping nan's
    if ~isempty(current_column(find(distances == min(distances), 1)))
        
        ideal_vector(1, current_item) = current_column(find(distances == min(distances), 1));
    
    else
        ideal_vector(1, current_item) = nan;
    end
    
end



%SMOOTH IDEAL_VECTOR
%
pre_smooth_vector(1,:) = smooth(ideal_vector(1,:), smooth_window);


%SHIFT IDEAL VECTOR BACK TO ORIGINAL min_max
%
%identify min_max corresponding to each padded row
min_maxs = nan(pads+1, 1);
min_maxs(:,1) = ones(pads+1,1).*min_max(1) + ones(pads+1,1).*((-pads/2)*min_max(2):min_max(2):(pads/2)*min_max(2))';

%fill ideal_vector(2,:) with the shift corresponding to which min_max the 
%corresponding item in ideal_vector(1,:) belongs to
for column = 1:length(pre_smooth_vector)
    pre_smooth_vector(2,column) = sum(pre_smooth_vector(1,column)>min_maxs) - ceil(size(padded_matrix, 1)/2);
end

%apply correction
smooth_vector = pre_smooth_vector(1,:) - (ones(1,length(pre_smooth_vector)).*(min_max(2).*pre_smooth_vector(2,:)));

end


