function [flags] = GetEV(filename, pos, starttime)
%GetEV(filename) takes a file path to a neurolynx .EV (flag/event) file as well
%as the position matrix (see GetVT.m) and outputs flags, a 6xN matrix.
%
% CURRENT OUTPUT VARIABLES:
%   event(:,1) = timestamps
%   event(:,2) = x coordinates
%   event(:,3) = y coordinates
%   event(:,4) = Event IDs
%   event(:,5) = TTLs
%   event(:,6) = Event Strings
%
% AVAILABLE OUTPUT VARIABLES:
%   TimeStamps: A 1xN vector of timestamps.
%   Event IDs: A 1xN vector of the event ID of each record.
%   TTLs: A 1xN vector of the decimal representation of the TTL value of each record.
%   Extras: A 8xN matrix of extra values for each record. These values are
%           generally not used.
%   Event Strings: A 1xN vector of the event string of each record.
%   Header: A Mx1 vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.


[TimeStamps, ~, ~, ~, ~, ~] = Nlx2MatEV_v3(filename, [1 1 1 1 1], 1, 1, []);

%places a 1 in rows with a timestamp (fifth column)
flags = [TimeStamps', nan(length(TimeStamps), 3), ones(length(TimeStamps), 1), NaN(length(TimeStamps), 1)];

if size(flags,1) > 2

    %converts timestamp units into seconds
    temp = ones(length(flags(:,1)), 1);
    starttimeEV(:,1) = temp.*starttime;
    flags(:,1) = flags(:,1) - starttimeEV(:,1);
    flags(:,1) = flags(:,1)./1000000;


    %add x and y coordinates to event by interpolating from position

    %interpx
    flags(:,2) = interp1(pos(:,1), pos(:,2), flags(:,1), 'linear');
    %interpy
    flags(:,3) = interp1(pos(:,1), pos(:,3), flags(:,1), 'linear');

    flags = flags(2:size(flags,1)-1, :);

else
    
    flags = [];
    
end

