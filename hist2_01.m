function matrix = hist2_01(eptrials, c, bins)
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
xedges = linspace(0, 1.0001, bins+1);
yedges = linspace(0, 1.0001, bins+1);


%filling xbin and ybin with firing rates. Last row is always 0, so we
%remove it. xn and yn are necessary. Don't ask me why.

%spikes
[xns, xbins] = histc(xs,xedges);
[yns, ybins] = histc(ys,yedges);
%time
[xnt, xbint] = histc(xt,xedges);
[ynt, ybint] = histc(yt,yedges);


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
matrix = matrix(1:end-1, 1:end-1);

%replace inf values with nans
matrix(matrix==inf)=nan;

%smooth
matrix = smooth2a(matrix, size(matrix,1)/10, size(matrix,2)/10);
matrix = smooth2a(matrix, size(matrix,1)/10, size(matrix,2)/10);

%fill in nans
matrix = inpaint_nans(matrix(1:bins, 1:bins));


        
%figure; 
%{ 
h = imagesc(1:length(matrix(:,1)),1:length(matrix(1,:)),matrix);
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
%}
