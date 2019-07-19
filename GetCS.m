function [CSC] = GetCS(filename, pos, starttime)

%GetCS(filename) takes a file path to a neurolynx .ncs (EEG) file and
%constructs a XXX
%
% Available outputs include: Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, and
% Header:
%
%   Notes on output data:
%   1. Each output variable's Nth element corresponds to the Nth element in
%      all the other output variables with the exception of the header output
%      variable.
%   2. The value of N in the output descriptions below is the total number of
%      records extracted.
%   3. For more information on Neuralynx records see:
%      http://www.neuralynx.com/static/software/NeuralynxDataFileFormats.pdf
%   4. Output data will always be assigned in the order indicated in the
%      FieldSelectionCSC. If data is not imported via a FieldSelectionCSC
%      index being 0, simply omit the output variable from the command.
%      EXAMPLE: FieldSelectionCSC = [1 0 0 0 1];
%      [Timestamps,Samples] = Nlx2MatCSC('test.ncs',FieldSelectionCSC,0,1,[]);
%
%   OUTPUT VARIABLES:
%   Timestamps: A 1xN vector of Timestamps.
%   ChannelNumbers: A 1xN vector of channel numbers.
%   SampleFrequencies: A 1xN vector of sample frequencies.
%   NumberOfValidSamples: A 1xN vector of the number of valid samples in the
%                         corresponding item in the Sample output variable.
%   Samples: A 512xN matrix of the data points. These values are in AD counts.
%   Header: A Mx1 vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.

[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC_v3(filename, [1 1 1 1 1], 1, 1, [] );

CSC_time = Timestamps';
CSC_samp = Samples;
CSC_samp = CSC_samp(:);

%converts timestamp units into seconds
temp = ones(length(CSC_time(:,1)), 1);
starttimeEV(:,1) = temp.*starttime;
CSC_time = CSC_time - starttimeEV(:,1);
CSC_time = CSC_time./1000000;
CSC_time(CSC_time<0, :) = [];

%"interpolate" timestamps
CSC_time = (CSC_time(1):((max(CSC_time)-CSC_time(1))/(length(CSC_samp)-512)):max(CSC_time))';
CSC_time = [CSC_time;zeros(511, 1)];

%restructure into larger matrix
CSC = [CSC_time, zeros(length(CSC_time), 2),  NaN(length(CSC_time), 1), zeros(length(CSC_time), 1), NaN(length(CSC_time), 1), CSC_samp];

%remove samples occuring after the last time stamp
CSC(find(CSC(:,1)==max(CSC(:,1))) + 1 :end, :) = [];

%DOWNSAMPLE to 2 kilohertz (1000 samples per second)
sampling_frequency = 1/(CSC(2,1)-CSC(1,1));
CSC = downsample(CSC,round(sampling_frequency/2000));

%add x and y coordinates to event by interpolating from pos
%interpx
CSC(:,2) = interp1(pos(:,1), pos(:,2), CSC(:,1), 'linear');
%interpy
CSC(:,3) = interp1(pos(:,1), pos(:,3), CSC(:,1), 'linear');

%Remove outer points to fit within the pos bounds?
CSC = CSC(2:length(CSC)-1, :);

end


