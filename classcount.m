function classcount(assigns_wb, bins, boi)
%makes a heatmap of spatial guesses in assigns_wb
%designed to work with ctxtspace_classify_crossbw
%
%boi is bin of interest

assigns_within = reshape(assigns_wb(:,1), 20, 36, 4);
assigns_between = reshape(assigns_wb(:,2), 20, 36, 4);



awthn=assigns_within(:, boi, :); 
awthn=awthn(:);
histc_within = histcounts(awthn, 1:bins^2+1)';
histc_within = histc_within./sum(histc_within(:));

figure; 
subplot(1,2,1)
imagesc(reshape(histc_within, bins, bins)); 
set(gca,'TickLength',[0, 0]); box off
colorbar; axis square
a = [0 .25];
caxis(a);
title within

[I_actual, J_actual] = ind2sub([bins, bins], find(reshape(1:bins^2, bins, bins)'==boi));
hold on; plot(I_actual, J_actual, 'ko', 'markersize', 30) %fix with ij coords of boi

abtwn=assigns_between(:, boi, :); 
abtwn=abtwn(:);
histc_between = histcounts(abtwn, 1:bins^2+1)';
histc_between = histc_between./sum(histc_between(:));


subplot(1,2,2) 
imagesc(reshape(histc_between, bins, bins));
set(gca,'TickLength',[0, 0]); box off
colorbar; axis square
caxis(a)
title between
hold on; plot(I_actual, J_actual, 'ko', 'markersize', 30) %fix with ij coords of boi



end