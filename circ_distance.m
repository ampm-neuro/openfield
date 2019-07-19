function [circ_dist, clockwise] = circ_distance(alpha, beta, min_max)
%circ_distance finds the circular distance between the two points"alpha" and
%"beta" that both exist within the circular range "range"
%
%e.g., 15 = circ_distance(5, 350, [0 360])
%
% clockwise is 1 if output distance is clockwise from alpha to beta, and -1
% distance is counter clockwise

%check for input error
if any(alpha < min_max(1)) || any(alpha > min_max(2)) || any(beta < min_max(1)) || any(beta > min_max(2))
    error('input values must be within input range')
end

%circular range
circ_range = abs(min_max(1) - min_max(2));

%non-circular distance between values
distance = abs(alpha - beta);

%find remainder after dividing by 360
remainder = rem(distance, repmat(circ_range,size(distance)));

%find circ_distance
circ_dist = nan(size(remainder));
idx = remainder > circ_range./2;
circ_dist(idx) = circ_range - remainder(idx);
circ_dist(~idx) = remainder(~idx);

%find direction
clockwise = nan(size(remainder));
clockwise((alpha-circ_dist)==beta) = -1;
clockwise(clockwise~=-1) = 1;

end