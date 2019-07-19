function [eptrials, clusters, context_ids, session_length] = loadnl(filename, varargin)
%function [pos, TimestampsVT, event, flags] = loadnl(filename)

%filename is the name of the folder (INSIDE THE MATLAB FOLDER) containing
%the day's data, e.g., 2013-10-02_08-05-23
%
%This outputs eptrials. See "trials_II" for more info.
%
%Requires that the data folder be inside a folder named 'neurodata' inside
%of the main MATLAB folder
%

%load excel files

context_ids = xlsread(strcat('/Volumes/ampm_PC2Mac/',num2str(filename),'/context_visits.xlsx'));
clusters = xlsread(strcat('/Volumes/ampm_PC2Mac/',num2str(filename),'/Cluster_descriptions.xlsx'));
if isequal(size(clusters), [2,1])
    clusters = [];
else
    clusters = clusters(~isnan(clusters(:,1)),:);
end




%Major neuralynx data files (video, spike, flag)


%USB KEY
%{
[pos, starttime, Targets, TimestampsVT, Angles] = GetVT(strcat('/Volumes/USB20FD/',num2str(filename),'/VT1.nvt'));
[event] = GetTT(strcat('/Volumes/USB20FD/',num2str(filename)), pos, starttime);
[flags] = GetEV(strcat('/Volumes/USB20FD/',num2str(filename),'/Events.nev'), pos, starttime);
%}

%MOBILE EXTERNAL HARD DRIVE
%
[pos, starttime] = GetVT(strcat('/Volumes/ampm_PC2Mac/',num2str(filename),'/VT1.nvt'));
event = GetTT(strcat('/Volumes/ampm_PC2Mac/',num2str(filename)), pos, starttime);
flags = GetEV(strcat('/Volumes/ampm_PC2Mac/',num2str(filename),'/Events.nev'), pos, starttime);
%}

%IMMOBILE EXTERNAL HARD DRIVE
%{
[pos, starttime, Targets, TimestampsVT, Angles] = GetVT(strcat('/Volumes/LaCie/PCAlte/',num2str(filename),'/VT1.nvt'));
[event] = GetTT(strcat('/Volumes/LaCie/PCAlte/',num2str(filename)), pos, starttime);
[flags] = GetEV(strcat('/Volumes/LaCie/PCAlte/',num2str(filename),'/Events.nev'), pos, starttime);
%}

%LOCAL DOCUMENTS
%{
[pos, starttime, TimestampsVT, angles, targets] = GetVT(strcat('/Users/ampm/Documents/MATLAB/PCia/Pcia_data/',num2str(subject), '/', num2str(filename),'/VT1.nvt'));
[event] = GetTT(strcat('/Users/ampm/Documents/MATLAB/PCia/Pcia_data/',num2str(subject), '/', num2str(filename)), pos, starttime);
[flags] = GetEV(strcat('/Users/ampm/Documents/MATLAB/PCia/Pcia_data/',num2str(subject), '/', num2str(filename),'/Events.nev'), pos, starttime);
%}

%THIS VARIABLE CONTROLS WHETHER THE CSC FILE IS ADDED TO eptrials
%if you would like the CSC data, it should be a 1. Otherwise, 0.
if nargin > 1
	INCLUDE_CSC = varargin{1};
else
    INCLUDE_CSC = 0;
end
if INCLUDE_CSC == 1
    CSC_hpc = GetCS(strcat('/Users/ampm/Documents/MATLAB/neurodata/',num2str(filename),'/CSC9.ncs'), pos, starttime);
    CSC_rsc = GetCS(strcat('/Users/ampm/Documents/MATLAB/neurodata/',num2str(filename),'/CSC10.ncs'), pos, starttime);
    
    eptrials = trials_pcia(event, pos, flags, context_ids, CSC);

elseif INCLUDE_CSC == 0

    eptrials = trials_pcia(event, pos, flags, context_ids);
    
end

session_length = max(eptrials(:,1))/60


if ~isempty(clusters)

    
    %check for common clusters mistake where multiple cells get the
    %same id
    if length(unique(eptrials(~isnan(eptrials(:,4)) & eptrials(:,4)>1, 4))) > length(clusters(:,1))
        previous_id = [];
        for ic = 1: size(clusters(:,1),1)
            if isequal(clusters(ic,1), previous_id)
                clusters(ic,1) = clusters(ic,1)+0.1;
            end
            previous_id = clusters(ic,1);
        end
    end 

    clusters(:,1) = clusters(:,1) + clusters(:,2)./10;
    clusters(:,2:3) = clusters(:,3:4); clusters(:,4) = [];

end

%}
end


