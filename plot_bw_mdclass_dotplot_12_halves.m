mkrsize = 30;
figure; 
colors = [11 178 27; 38 166 230; 202 16 16; 228 146 40]./255;
close

%FIRST TWO AXES
figure; hold on; axis([.65 1.67 .65 1.67]); axis square; set(gca,'TickLength',[0, 0]); box off
xlabel('distance to black')
ylabel('distance to white')
%each time window
axes = [1 2];
xplot = []; yplot = [];
for vis = unique(ctx_context_bw)'
    %hold
    xplot = [xplot distances_bw(class_test_bw == vis, axes(1))];
    yplot = [yplot distances_bw(class_test_bw == vis, axes(2))];
end



for tw = 1:size(xplot,1)
   for vis = 1:size(xplot,2)
      plot(xplot(tw, vis), yplot(tw,vis), '.', 'markersize', mkrsize, 'linewidth', 1, 'color', colors(vis,:)) 
   end
end

%
for vis = 1:size(xplot,2)
        %mean
        plot(mean(xplot(:, vis)), mean(yplot(:,vis)), '.', 'markersize', 3*mkrsize, 'linewidth', 5, 'color', colors(vis,:))
        plot(mean(xplot(:, vis)), mean(yplot(:,vis)), 'o', 'markersize', 3*mkrsize/pi, 'linewidth', 5, 'color', [0 0 0])
end
%}

plot([0 2], [0 2], 'k--', 'linewidth', 4)
hold on; plot([1 1], [0 2], 'k--', 'linewidth', 1)
hold on; plot([0 2], [1 1], 'k--', 'linewidth', 1)