function avg_circ_vect = circ_mean(input_vector, min_max)

%make sure input_vector is a row vector
if size(input_vector,1) > size(input_vector,2)
    input_vector = input_vector';
end

%SOME PRE-CODING THINGS
%

%pad_vector is a local function that calculates the pad vector by shifting 
%input_vector up or down
function shifted_vector = shift_vector(vector, shift, min_max)
    
    min_max_shift = ones(1,length(vector)).*(min_max(2)*shift);
    shifted_vector = vector + min_max_shift;

end



%BUILD PADDED_MATRIX
%

%pads determines the number of rows that surround input_vector to make 
%the padded_matrix. Must be an even number.
pads = 100;

%preallocate padded_matrix
padded_matrix = nan(pads+1, length(input_vector));

%build padded_matrix with shifted vectors
for row = 1 : size(padded_matrix, 1)
    
    %pad_vector is a local function
    padded_matrix(row,:) = shift_vector(input_vector, row - ceil(size(padded_matrix, 1)/2), min_max);

end



%CALCULATE IDEAL_VECTOR
%

%preallocate ideal_vector
ideal_vector = nan(1, length(input_vector));

%start with the first item in input_vector
ideal_vector(1,1) = input_vector(1);

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



%TAKE MEAN OF IDEAL_VECTOR
%
ideal_average = mean(ideal_vector(1,:));


%SHIFT IDEAL_AVERAGE back to within bounds defined by min_max
%

if ideal_average > min_max(1) && ideal_average < min_max(2)
    
    avg_circ_vect = ideal_average;
    
elseif ideal_average < min_max(1)
    
    avg_circ_vect = min_max(2) + rem(ideal_average, min_max(2));
    
elseif ideal_average > min_max(2)
    
  	avg_circ_vect = rem(ideal_average, min_max(2));

    
end


avg_circ_vect(avg_circ_vect>min_max(2)) = avg_circ_vect(avg_circ_vect>min_max(2))-min_max(2);

end







