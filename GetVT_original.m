function [pos, starttime, TimestampsVT, angles, targets] = GetVT(filename)
%GetVT(filename) takes a file path to a neurolynx .nvt (video) file and
%constructs a 4xN position matrix comprised of columns: Timestamps, X, Y, and
%head direction. The X and Y vectors are then re-written with interpolated
%data (to fill in missed elements). The position matrix, pos, is then output.
%
% Available vectors include: Timestamp, X, Y, Angle, Targets, Points, and
% Header:
%   TimestampsVT: A 1xN vector of timestamps.
%   Extracted X: A 1xN vector of the calculated X coordinate for each record.
%   Extracted Y: A 1xN vector of the calculated Y coordinate for each record.
%   Extracted Angle: A 1xN vector of the calculated head direction angle for
%                    each record. This value is in degrees.
%   Targets: A 50xN matrix of the targets found for each frame. These values
%            are encoded using the VT bitfield encoding.
%   Points: A 480xN matrix of the threshold crossings found for each frame.
%           These values are encoded using the VT bitfield encoding.
%   Header: A Mx1 vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.


%ESTABLISH pos FROM NEURALYNX OUTPUT
%

%output from Nlx2MatVT_v3, the mac version of standard neuralynx extraction
%code

%filename
[TimestampsVT, X, Y, angles, targets, ~, ~] = Nlx2MatVT_v3(filename, [1 1 1 1 1 1], 1, 1, []);

%reverse Y, which fixes a camera thing
Y = abs(Y-ones(size(Y))*max(Y));

%translate outputs into pos, the begining of the position vector that is
%input into trials_III
pos = [TimestampsVT', X', Y', NaN(length(TimestampsVT), 1), zeros(length(TimestampsVT), 1), angles'];


%INTERPOLATE X and Y vectors
%

%builds vector of timestamps that have associated x and y values
tfull = pos(pos(:,2)>0 & pos(:,3)>0,1);

%builds vectors of the x and y values associated with above time stamps
pfullx = pos(pos(:,2)>0 & pos(:,3)>0,2);
pfully = pos(pos(:,2)>0 & pos(:,3)>0,3);

%angle vector associated with above time stamps
afull = pos(pos(:,2)>0 & pos(:,3)>0,6);
afull(afull==0) = nan;

pos = [tfull pfullx pfully ones(size(tfull)) zeros(size(tfull)) afull];

%converts timestamp units into seconds
%{
starttime = zeros(length(tfull(:,1)), 1);
starttime(:,1) = tfull(1,1);
tfull(:,1) = tfull(:,1) - starttime(:,1);
tfull(:,1) = tfull(:,1)./1000000;


%resetting pos. FIX THIS IF TRYING TO INCORPORATE MORE ITEMS FROM NEURALYNX
pos = [NaN(floor(max(tfull(:,1)))*100 + 1, 3), ones(floor(max(tfull(:,1)))*100 + 1, 1), zeros(floor(max(tfull(:,1)))*100 + 1, 1), NaN(floor(max(tfull(:,1)))*100 + 1, 1)];

%generating desired timestamps (10ms increments)
pos(:,1) = (0:0.01:floor(max(tfull(:,1))))';

%re-writes the position X and Y vectors with interpolated data
pos(:,2) = interp1(tfull ,pfullx, pos(:,1), 'linear');
pos(:,3) = interp1(tfull ,pfully, pos(:,1), 'linear');


%re-wrires the HD angles to match up with new positions
pos(pos(:,6)==0,6) = nan; %delete dropped signals
pos(:,6) = circ_interp1(afull(~isnan(afull))', [0 360], pos(:,1)', tfull(~isnan(afull))')';

%starttime
starttime = starttime(1,1);
%}






