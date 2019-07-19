function [interped_vector] = circ_interp1(subsampled_vector, min_max, complete_vector, partial_vector)
%
% function [interp_vector] = circ_interp1(subsampled_vector, min_max, complete_vector, partial_vector)
%
% circ_interp1 applies the matlab interp1 function on circular data, and is
% closely based on the code circ_smooth by ampm
%
% interp1 (see help interp1) finds interped_vector, which is subsampled_vector 
% with new values at all the points in complete_vector that are not also in
% partial_vector. Essentially it finds interped_vector where:
%
% subsampled_vector : interped_vector  ::  partial_vector : complete_vector
%
% e.g., head_direction points you have : all the head direction points
% :: time points at which you have head direction points : all the
% time points
%
% Dealing with circularity:
% circ_interp1 forgoes any knowledge about circular mathematics in favor of 
% brute force. Hopefully you have a shiny computer.
%
% The general strategy is to pad subsampled_vector above the max and below the min 
% with copies of itself. An ideal vector is then calculated by iterating through 
% the columns of the padded_matrix and selecting the item in the current column 
% with the minimum distance to the selected item in the previous column. The interp1
% function is then run on the ideal vector. We finish by returning the items 
% in the ideal vector back to their original position within the range defined 
% by min_max.
%
% Inputs:
%   ALL INPUTS ARE ROW VECTORS
%
%  subsampled_vector is the vector of items that you want to expand with
%   interplated points
%
%  min_max describes the way in which the data in subsampled_vector is
%   circular by giving the minumum and maximum values of the circle such that
%   the minimum value is equal to or at contiguous with the maximum value 
%   (e.g., for degrees around a circle, min_max = [0 360]).
%
%  complete_vector is the vector containing all the items at which you
%   would like to expand subsampled_vector
%
%  partial_vector is a subset of complete_vector that only contains the
%   items corresponding to values that already exist in subsampled_vector
%
% Output:
%  interped_vector is the interpolated vector and contains all the items in 
%   subsampled_vector as well a number of new points that is equal to the
%   length(complete_vector) - length(partial_vector)
%

%pad_vector is a local function that calculates the pad vector by shifting 
%subsampled_vector up or down
function shifted_vector = shift_vector(vector, shift, min_max)
    min_max_shift = ones(1,length(vector)).*(min_max(2)*shift);
    shifted_vector = vector + min_max_shift;
end


%BUILD PADDED_MATRIX
%

%pads determines the number of rows that surround subsampled_vector to make 
%the padded_matrix. Must be an even number.
pads = 1000;

%preallocate padded_matrix
padded_matrix = nan(pads+1, length(subsampled_vector));

%build padded_matrix with shifted vectors
for row = 1 : size(padded_matrix, 1)
    
    %row
    %ceil(size(padded_matrix, 1)/2)
    
    %pad_vector is a local function
    
    padded_matrix(row,:) = shift_vector(subsampled_vector, row - ceil(size(padded_matrix, 1)/2), min_max);

end
%CALCULATE IDEAL_VECTOR
%

%preallocate ideal_vector
ideal_vector = nan(1, length(subsampled_vector));

%start with the first item in subsampled_vector
ideal_vector(1,1) = subsampled_vector(1);

%find the next item
for current_item = 2:length(ideal_vector)
        
    %indexing out just the column under consideration
    current_column = padded_matrix(:, current_item);
    
    %a column full of the current item
    temp_item_column = ones(size(current_column)).*ideal_vector(1, current_item-1);
    
    %how far is each item in the current column from the current item?
    distances = abs(current_column - temp_item_column);
    
    %fill ideal_vector with the closest item and then shift the corresponds 
    %to the row it came from
    ideal_vector(1, current_item) = current_column(find(distances == min(distances), 1));
    
end

%Interp1 IDEAL_VECTOR
%

%establish ingredients for interp1
partial_x = partial_vector; %time_vt
partial_y = ideal_vector(1,:);
complete_x = complete_vector; %eptrials(eptrials(:,14) == 1,1)'

%see help interp1, but linear is the most straightforward
method = 'linear';

%unique items only
[partial_x, idx_unique_partial_x] = unique(partial_x);
partial_y = partial_y(idx_unique_partial_x);

%interp1
complete_y = interp1(partial_x, partial_y, complete_x, method); %interp1

%ouput interpolated ideal_vector
ideal_vector = complete_y;

%SHIFT IDEAL VECTOR BACK TO ORIGINAL min_max
%

%identify min_max corresponding to each padded row
%min_maxs = nan(pads+1, 1);
%min_maxs(:,1) = ones(pads+1,1).*min_max(1) + ones(pads+1,1).*((-pads/2)*min_max(2):min_max(2):(pads/2)*min_max(2))';

%fill ideal_vector(2,:) with the shift corresponding to which min_max the 
%corresponding item in ideal_vector(1,:) belongs to. (This is on the clever side)
%for column = 1:length(ideal_vector)
%    ideal_vector(column) = sum(ideal_vector(1,column)>min_maxs) - ceil(size(padded_matrix, 1)/2);
%end

ideal_vector(ideal_vector>0) = rem(ideal_vector(ideal_vector>0),360); 
ideal_vector(ideal_vector<0) = rem(ideal_vector(ideal_vector<0),360) + 360; 

%apply correction
interped_vector = ideal_vector;%(1,:) - (ones(1,length(ideal_vector)).*(min_max(2).*ideal_vector(2,:)));

end

