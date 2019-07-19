function visit_matrix = position_heatmap(eptrials, bins)
% function position_heatmap(eptrials, c, bins)
%
% Extract 2D histogram data containing the positions each
% of the bins defined by the x and y ranges in eptrials, and the bin
% spacing defined by bins.


%samplingrate
smplrt=length(eptrials(eptrials(:,4)==1,1))/(max(eptrials(:,1)) - min(eptrials(:,1)));

%all time samples in two vectors 
xt = eptrials(eptrials(:,4)==1, 2);
yt = eptrials(eptrials(:,4)==1, 3);

%evenly spaced bins of x and y coordinate ranges (incl pos - not just event -
%data)
xedges = linspace(min(eptrials(:,2)), max(eptrials(:,2)), bins);
yedges = linspace(min(eptrials(:,3)), max(eptrials(:,3)), bins);

%majik
[~, xbint] = histc(xt,xedges);
[~, ybint] = histc(yt,yedges);

xnbint = length(xedges);
ynbint = length(yedges);
xbint(xbint == 0) = 1;
ybint(ybint == 0) = 1;

if xnbint >= ynbint
    xyt = ybint*(xnbint) + xbint;
    indexshiftt = xnbint;
else
    xyt = xbint*(ynbint) + ybint;
    indexshiftt = ynbint;
end

xyunit = unique(xyt);
hstrest = histc(xyt,xyunit);

%establish the histmat matrix
histmatt = zeros(xnbint, ynbint);
histmatt(xyunit-indexshiftt) = hstrest;
histmatt = histmatt'./smplrt;
histmatt(1,1)=0;

%smooth
size_fix = histmatt; 
size_fix(~isnan(size_fix))=1;
mask = [1 3 1; 3 4 3; 1 3 1]./20;
histmatt = conv2nan(histmatt, mask, 'same');histmatt = histmatt.*size_fix;

%output
visit_matrix = histmatt./sum(histmatt(:));

%plots output from hist2
%{
figure; 
pcolor(1:length(visit_matrix(:,1)),1:length(visit_matrix(1,:)),visit_matrix);
colorbar; 
axis square tight;
set(gca,'ydir','reverse')
shading flat
%}

end