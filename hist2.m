function matrix = hist2(eptrials, c, bins)
% function histmat  = hist2(eptrials, c, bins)
%
% Extract 2D histogram data containing the firing rate for cell c at each
% of the bins defined by the x and y ranges in eptrials, and the bin
% spacing defined by bins.
%
% Essentially this builds two 2d histograms in parrallel, and then plots one
% divided by the other.
%
% This is a modified version of Dave Bulkin's hist2.m.

%samplingrate
smplrt=length(eptrials(eptrials(:,4)==1,1))/(max(eptrials(:,1)) - min(eptrials(:,1)));

%subset of trials
%eptrials = eptrials(eptrials(:,5)>20, :);

%all spike events in two vectors
xs = eptrials(eptrials(:,4)==c, 2);
ys = eptrials(eptrials(:,4)==c, 3);

%all time samples in two vectors 
xt = eptrials(eptrials(:,4)==1, 2);
yt = eptrials(eptrials(:,4)==1, 3);

%evenly spaced bins of x and y coordinate ranges (incl pos - not just event -
%data)
xedges = linspace(min(eptrials(:,2)), max(eptrials(:,2)), bins);
yedges = linspace(min(eptrials(:,3)), max(eptrials(:,3)), bins);


%filling xbin and ybin with firing rates. Last row is always 0, so we
%remove it. xn and yn are necessary. Don't ask me why.

%spikes
[xns, xbins] = histc(xs,xedges);
[yns, ybins] = histc(ys,yedges);
%time
[xnt, xbint] = histc(xt,xedges);
[ynt, ybint] = histc(yt,yedges);


%THIS SECTION REMOVES PIXLES THAT WERE ONLY VISITED ON ONE TRIAL. IT IS NOT
%APPROPRIATE FOR NON-TRIAL BASED DATA (CAN BE COMMENTED OUT).
%{
    %matrix y-cords x-cords and trial number
    pixles_and_trials = [ybint xbint eptrials(eptrials(:,10)==1 & eptrials(:,8)==1, 5)]; 

    %get indices of unique rows (each pixle visited on each trial)
    [~, uni_indi, ~] = unique(pixles_and_trials, 'rows');
    pixle_per_trial = pixles_and_trials(uni_indi, 1:2);

    %get repeated rows of pixles_and_trials(indices, 1:2) (pixles that were 
    %visited on multiple trials)
    [~,~,n] = unique(pixle_per_trial, 'rows'); %see help unique
    hist_counts = hist(n, .5:1:(max(n)-.5)); %how many trials each pixle was visited on
    multi_trial_pixles = pixle_per_trial(ismember(n, find(hist_counts>1)), :); %pixles visited on >1 trials

    %index pixles_and_trials for pixles visited on multiple trials
    multivisit_pixle_index_space = logical(ismember([ybint xbint], multi_trial_pixles, 'rows'));
    multivisit_pixle_index_event = logical(ismember([ybins xbins], multi_trial_pixles, 'rows'));

    %re-define set of visited x and y matrix coords to only include
    %pixles visited on multiple trials
    xbint = xbint(multivisit_pixle_index_space);
    ybint = ybint(multivisit_pixle_index_space);
    xbins = xbins(multivisit_pixle_index_event);
    ybins = ybins(multivisit_pixle_index_event);
%}

%xbin, ybin zero for out of range values (see the help of histc) force this 
%event to the first bins

%spikes
xbins(find(xbins == 0)) = 1;
ybins(find(ybins == 0)) = 1;
xbint(find(xbint == 0)) = 1;
ybint(find(ybint == 0)) = 1;
%time
xnbins = length(xedges);
ynbins = length(yedges);
xnbint = length(xedges);
ynbint = length(yedges);

%spikes
if xnbins >= ynbins
    xys = ybins*(xnbins) + xbins;
    indexshifts = xnbins;
else
    xys = xbins*(ynbins) + ybins;
    indexshifts = ynbins;
end

%time
if xnbint >= ynbint
    xyt = ybint*(xnbint) + xbint;
    indexshiftt = xnbint;
else
    xyt = xbint*(ynbint) + ybint;
    indexshiftt = ynbint;
end


%spikes
xyunis = unique(xys);
hstress = histc(xys,xyunis);
%time
xyunit = unique(xyt);
hstrest = histc(xyt,xyunit);

%establish the histmat matrix
%spikes
histmats = zeros(xnbins, ynbins);
%time
histmatt = zeros(xnbint, ynbint);

%spikes
histmats(xyunis-indexshifts) = hstress;
histmats = histmats';
%time
histmatt(xyunit-indexshiftt) = hstrest;
histmatt = histmatt'./smplrt;

%spikes
histmats(1,1)=0;
%time
histmatt(1,1)=0;

%plot
matrix = histmats./histmatt;


%replace inf values
if sum(isinf(matrix(:)))>0
%matrix(isinf(matrix))=NaN;%replace inf values with NaN

    %set inf pixles to nanmean of surrounding pixles
    [inf_y, inf_x] = ind2sub(size(matrix), find(matrix==inf));

        %remove corners (if they exist)
        try
            inf_corners = [1 1; size(matrix, 1) 1; 1 size(matrix, 2); size(matrix)];
            inf_y = inf_y(~ismember([inf_y inf_x], inf_corners, 'rows'));
            inf_x = inf_x(~ismember([inf_y inf_x], inf_corners, 'rows'));
        catch 
            inf_y 
            inf_x
        end
    
    for inf_value = 1:length(inf_y)
        
        %gather surrounding pixle firing rates
        if inf_y(inf_value) == 1
         	surround = [matrix(inf_y(inf_value)+1, inf_x(inf_value)) matrix(inf_y(inf_value), inf_x(inf_value)+1) matrix(inf_y(inf_value), inf_x(inf_value)-1) matrix(inf_y(inf_value)+1, inf_x(inf_value)+1) matrix(inf_y(inf_value)+1, inf_x(inf_value)-1)];
        elseif inf_y(inf_value) == size(matrix,1)
         	surround = [matrix(inf_y(inf_value)-1, inf_x(inf_value)) matrix(inf_y(inf_value), inf_x(inf_value)+1) matrix(inf_y(inf_value), inf_x(inf_value)-1) matrix(inf_y(inf_value)-1, inf_x(inf_value)-1) matrix(inf_y(inf_value)-1, inf_x(inf_value)+1)];
        elseif inf_x(inf_value) == 1
         	surround = [matrix(inf_y(inf_value)+1, inf_x(inf_value)) matrix(inf_y(inf_value)-1, inf_x(inf_value)) matrix(inf_y(inf_value), inf_x(inf_value)+1) matrix(inf_y(inf_value)+1, inf_x(inf_value)+1) matrix(inf_y(inf_value)-1, inf_x(inf_value)+1)];
        elseif inf_x(inf_value) == size(matrix,2)
          	surround = [matrix(inf_y(inf_value)+1, inf_x(inf_value)) matrix(inf_y(inf_value)-1, inf_x(inf_value)) matrix(inf_y(inf_value), inf_x(inf_value)-1) matrix(inf_y(inf_value)-1, inf_x(inf_value)-1) matrix(inf_y(inf_value)+1, inf_x(inf_value)-1)];
        else
            surround = [matrix(inf_y(inf_value)+1, inf_x(inf_value)) matrix(inf_y(inf_value)-1, inf_x(inf_value)) matrix(inf_y(inf_value), inf_x(inf_value)+1) matrix(inf_y(inf_value), inf_x(inf_value)-1) matrix(inf_y(inf_value)+1, inf_x(inf_value)+1) matrix(inf_y(inf_value)-1, inf_x(inf_value)-1) matrix(inf_y(inf_value)+1, inf_x(inf_value)-1) matrix(inf_y(inf_value)-1, inf_x(inf_value)+1)];
        end
        
        if ~isempty(surround)
        
            %set new value
            matrix(inf_y(inf_value), inf_x(inf_value)) = nanmean(surround(~isinf(surround)));
        
        else
        
            matrix(inf_y(inf_value), inf_x(inf_value)) = NaN;
            
        end
        
    end
    
   %on the very unlikely case that we actually removed corners - replace
   %with mean overall rate
   matrix(isinf(matrix))=mean(matrix(~isinf(matrix)));%replace remaining inf values with mean firing rate
   
end


%matrix with NaNs in place of NaNs and 1's in place of numbers. This
%prevents (actually, deletes) the bloat occuring during convolution.
size_fix = matrix; 
size_fix(~isnan(size_fix))=1;


mask = [1 3 1; 3 4 3; 1 3 1]./20;
matrix = conv2nan(matrix, mask, 'same');%matrix = matrix.*size_fix;
%matrix = conv2nan(matrix, mask, 'same');%matrix = matrix.*size_fix;
%matrix = conv2nan(matrix, mask, 'same');%matrix = matrix.*size_fix;
%matrix = conv2nan(matrix, mask, 'same');%matrix = matrix.*size_fix;
matrix = conv2nan(matrix, mask, 'same');matrix = matrix.*size_fix;
%matrix = conv2nan(matrix, mask, 'same');matrix = matrix.*size_fix;

%matrix = interp2(matrix, 3);

%krnl=gausswin(5)*gausswin(5)'; % create a gausswin kernel
%mask=krnl./sum(sum(krnl)); %normalize the kernel, so it doesn't alter firing rate
%matrix = conv2(matrix, mask, 'same');

%plots output from hist2

%figure; 
%pcolor(xedges,yedges,matrix); 
h = pcolor(1:length(matrix(:,1)),1:length(matrix(1,:)),matrix);
colorbar; 
axis square tight;
set(gca,'ydir','reverse')
shading flat

%higher numbers values for 'a' will result in warmer colors overall
a = 1.1;
b = 1.1;

clrrng(1) = min(matrix(:));
clrrng(2) = max(matrix(:));

e = clrrng(1)*b;
d = clrrng(2)/a;

caxis([e d])
title(['Cell ',num2str(c)],'fontsize', 16) 
