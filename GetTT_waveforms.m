function [correlations, FRs_Amp, thresholds, peaks, all_cell_idx, amplitudes] = GetTT_waveforms(filename, cluster_idx, session_file)
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

    if exist(strcat(num2str(filename), '\TT', num2str(i), '_sorted.NTT')) > 0
        
        a=a+1;
        
        b{a} = strcat(num2str(filename), '\TT', num2str(i), '_sorted.NTT');
        
        c(a) = i;
    end

end

%fills cell array 'events' with the matrices extracted from each sorted .NVT file.
all_cell_idx = [];
thresholds = [];
correlations = [];
FRs_Amp = [];
amplitudes = cell(a, 1);

%cluster index counter
cic = 0;


%if there are clusters
if a>0
    
    events = cell(length(b),1);

    %each tt
    for i = 1:length(b)
   
        %extract data
        %[TimestampsTT, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike_v3(num2str(b{i}), [1 1 1 1 1], 1, 1, [] );
        [TimestampsTT,CellNumbers, Features] = Nlx2MatSpike(num2str(b{i}), [1 0 1 1 0], 0, 1, [] );

        %Features_orig = Features;
        %threshold
        threshold_hold = min(max(Features));
        
        %fix time
        TimestampsTT = (TimestampsTT - min(TimestampsTT))./100000;
        
        %remove everything associated with unsorted events
        TimestampsTT(CellNumbers==0)=[];
        Features(:, CellNumbers==0)=[];
        CellNumbers(CellNumbers==0)=[];
        
        for each_cell = unique(CellNumbers)
            cic = cic+1;
            if cluster_idx(cic) == 0
                continue
            end

            %calculate firing rates (30s bins)
            [FRs, edges] = histcounts(TimestampsTT(CellNumbers==each_cell), floor(max(TimestampsTT(CellNumbers==each_cell))/30));

            %restrict Features to peaks (amplitude)
            peaks = max(Features); peaks = peaks(CellNumbers==each_cell);
            mean_amplitudes = nan(1, length(edges)-1);
            for itw = 1:length(edges)-1
                c_edges = edges([itw itw+1]);
                mean_amplitudes(itw) = mean(peaks(TimestampsTT(CellNumbers==each_cell) >= c_edges(1) & TimestampsTT(CellNumbers==each_cell) < c_edges(2)));
            end

            %calculate correlation outputs
            [pearsonr, pval] = corr(FRs(~isnan(mean_amplitudes))', mean_amplitudes(~isnan(mean_amplitudes))');
            
            %fit_line(FRs',mean_amplitudes');
            thresholds = [thresholds; [threshold_hold min(peaks)]];
            correlations = [correlations; [pearsonr, pval]];

            FRA_hold = [FRs; mean_amplitudes];
            FRs_Amp = [FRs_Amp;  mat2cell(FRA_hold, size(FRA_hold, 1), size(FRA_hold, 2))];
            amplitudes{i} = peaks;
            
            %plot!
            %{
                plot_peak_amp(peaks, threshold_hold, 1:length(peaks)); set(gcf,'pos',[0 600 900 500])

                %error yo

                choice = questdlg('Is this a good cell?','Plot Check','Yes','No','No');
                switch choice
                    case 'Yes'
                        all_cell_idx = [all_cell_idx; 1];
                        print(strcat('/Users/ampm/Desktop/papers/Context/pcia_figs/d5_figs/time_sort_select/', num2str(session_file),'_tt_' ,num2str(i),'_cell_' ,num2str(each_cell), '.pdf'), '-dpdf')
                        close
                    otherwise
                        all_cell_idx = [all_cell_idx; 0];
                        close
                end
            %}
        
        end
    end


else
    
    correlations = [];
    warning('No clusters identified')

end



