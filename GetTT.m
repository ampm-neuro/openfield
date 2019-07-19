function [event] = GetTT(filename, pos, starttime)
%GetTT(filename) takes a file path to a neurolynx .NTT (event) file as well
%as the position matrix (see GetVT.m) and outputs event, a 4xN matrix.
%
% CURRENT OUTPUT VARIABLES:
%   event(:,1) = timestamps
%   event(:,2) = x coordinates
%   event(:,3) = y coordinates
%   event(:,4) = cluster ids
%   event(:,5) = tetrode of origin
%
% AVAILABLE OUTPUT VARIABLES:
%   TimestampsTT: A 1xN vector of timestamps.
%   ScNumbers: A 1xN vector of spike channel numbers. This is the order that
%              the spike AEs were created and have nothing to do with the
%              AD channel number.
%   CellNumbers: A 1xN vector of classified cell numbers. If no cell was
%                classified for this spike, this value will be zero.
%   Features: A 8xN vector of the features (e.g. Peak, Valley, etc.) calculated
%             by Cheetah.
%   Samples: A 32xMxN matrix of the data points. Where M is the number of
%            subchannels in the spike file (NTT M = 4, NST M = 2, NSE M = 1).
%            These values are in AD counts.
%   Header: A Mx1 vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.




%identifies sorted .NTT files and orders their file paths in 'b', a cell
%array. c keeps track of the TT numbers.

a=0;

for i = 1:16

    if exist(strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT')) > 0 %|| exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
        
        a=a+1;
        
        %if exist(strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT')) > 0
        
            b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT');
        
        %elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
            
            %b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT');
            
        %end
        
        c(a) = i;
    
   % elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_SS_01.NTT')) > 0
        
    elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedanna.NTT')) > 0 %|| exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
        
        a=a+1;
        
        %if exist(strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT')) > 0
        
            b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sortedanna.NTT');
        
        %elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
            
            %b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT');
            
        %end
        
        c(a) = i;    
        
    elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0 %|| exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
        
        a=a+1;
        
        %if exist(strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT')) > 0
        
            b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT');
        
        %elseif exist(strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT')) > 0
            
            %b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sortedAnna.NTT');
            
        %end
        
        c(a) = i;    
        
    end

end

%fills cell array 'events' with the matrices extracted from each sorted .NVT file.

%if there are clusters
if a>0
    events = cell(length(b),1);

    for i = 1:length(b)
   
        [TimestampsTT, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike_v3(num2str(b{i}), [1 1 1 1 1], 1, 1, [] );

        d = [TimestampsTT', zeros(length(TimestampsTT), 2), CellNumbers', zeros(length(TimestampsTT), 1), NaN(length(TimestampsTT), 1)];
    
        %removes timestamps (entire rows) associated with noise (cluster 0)
        d(d(:,4)==0, :)=[];
    
        %lables each cell as TT_of_orgin.cluster_ID (eg TT 1 cluster 1 is 1.0100)
        d(:,4) = d(:,4)./10 + repmat(c(i), size(d(:,4)));

        %fills cell array
        events{i} = d;
    
    end

%merges the many matrices into a single matrix sorted by timestamps
event = sortrows(cell2mat(events), 1);

%converts timestamp units into seconds
temp = ones(length(event(:,1)), 1);
starttimeTT(:,1) = temp.*starttime;
event(:,1) = event(:,1) - starttimeTT(:,1);
event(:,1) = event(:,1)./1000000;

%adds x and y coordinates to event by interpolating from position
%interpx
event(:,2) = interp1(pos(:,1), pos(:,2), event(:,1), 'linear');
%interpy
event(:,3) = interp1(pos(:,1), pos(:,3), event(:,1), 'linear');

else
    
    event = [];
    warning('No clusters identified')
    
end



