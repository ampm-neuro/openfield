function [FRs] = GetTT_waveforms2(filename, cluster_idx, session_file, aci)
%GetTT(filename) takes a file path to a neurolynx .NTT (event) file as well
%as the position matrix (see GetVT.m) and outputs event, a 4xN matrix.
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

    if exist(strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT')) > 0
        
        a=a+1;
        
        b{a} = strcat(num2str(filename), '/TT', num2str(i), '_sorted.NTT');
        
        c(a) = i;
    end

end

%fills cell array 'events' with the matrices extracted from each sorted .NVT file.
FRs = [];


%cluster index counter
cic = 0;
cic2 = 0;


%if there are clusters
if a>0
    
    events = cell(length(b),1);

    %each tt
    for i = 1:length(b)
   
        %extract data
        %[TimestampsTT, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike_v3(num2str(b{i}), [1 1 1 1 1], 1, 1, [] );
        [TimestampsTT,CellNumbers, Features] = Nlx2MatSpike_v3(num2str(b{i}), [1 0 1 1 0], 0, 1, [] );
        [TimeStampsEV, ~, ~, ~, ~, ~] = Nlx2MatEV_v3([filename '/Events.nev'], [1 1 1 1 1], 1, 1, []);
        %Features_orig = Features;
        %threshold
        threshold_hold = min(max(Features));
        
        %fix time
        TimestampsTT = (TimestampsTT - min(TimestampsTT))./1000000;
        TimeStampsEV = (TimeStampsEV - min(TimeStampsEV))./1000000;
        %remove everything associated with unsorted events
        TimestampsTT(CellNumbers==0)=[];
        Features(:, CellNumbers==0)=[];
        CellNumbers(CellNumbers==0)=[];
        
        for each_cell = unique(CellNumbers)
            %check if meets cluster constraints
            cic = cic+1;
            if cluster_idx(cic) == 0
                continue
            end
            %check if previously selected good amplitude
            cic2 = cic2+1;
            if aci(cic2) == 0
                continue
            end
            
            %calculate firing rates
            binsize = 2;%seconds
            [FRs_local] = histcounts(TimestampsTT(CellNumbers==each_cell), floor(max(TimestampsTT(CellNumbers==each_cell))/binsize));
            FRs_local = FRs_local./binsize; %set to hz
            
            %restrict Features to peaks (amplitude)
            peaks = max(Features); peaks = peaks(CellNumbers==each_cell);
            peaktimes = TimestampsTT(CellNumbers==each_cell);
            %load
            FRs = [FRs;  {FRs_local}];
            
            
            
            %{
            %plot amplitudes
            figure;
            plot_peak_amp(peaks, threshold_hold, peaktimes); 
            set(gcf,'pos',[0 600 900 500])
            
            title([num2str(session_file),'-tt-' ,num2str(i),'-cell-' ,num2str(each_cell)])
            
            %overlay events
            hold on
            for ev = 1:length(TimeStampsEV)
                plot([TimeStampsEV(ev) TimeStampsEV(ev)], ylim, 'r-')
            end
            hold off
            xticklabels(xticks/60)
            set(gca,'TickLength',[0, 0]); box off
            
            %print(strcat('/Users/ampm/Desktop/papers/Context/pcia_figs/d5_figs/time_sort_select/',...
            %    num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell), 'amp.pdf'), '-dpdf', '-painters')

            print(strcat('/Users/ampm/Desktop/papers/Context/pcia_figs/d5_figs/all_rate_examples/',...
                num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell), 'amp.pdf'),'-bestfit', '-dpdf')
            %}
            
            %plot rates
            figure; 
            fit_line_ratetime(FRs_local, binsize);
            set(gcf,'pos',[0 600 900 500])
            set(gca,'TickLength',[0, 0]); box off
            
            title([num2str(session_file),'-tt-' ,num2str(i),'-cell-' ,num2str(each_cell)])
            xlabel('Time from Start (min)')
            ylabel('Firing Rate (Hz)')
            
            %events
            hold on
            for ev = 1:length(TimeStampsEV)
                plot([TimeStampsEV(ev) TimeStampsEV(ev)], ylim, 'r-')
            end
            hold off
            xticklabels((xticks)/60)
            
            %print(strcat('/Users/ampm/Desktop/papers/Context/pcia_figs/d5_figs/time_sort_select/',...
            %    num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell), 'rate.pdf'), '-dpdf', '-painters')
            
            print(strcat('/Users/ampm/Desktop/papers/Context/pcia_figs/d5_figs/all_rate_examples/',...
                num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell), 'rate.pdf'),'-bestfit', '-dpdf')
            %}
        
            close all
            
        end
    end


else
    
    correlations = [];
    warning('No clusters identified')

end



