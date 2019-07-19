function plot_bw_mdclass(distances_bw_1, distances_bw_2, class_test_1, class_test_2) 

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

%first half
xplot_1 = []; yplot_1 = [];
for vis = unique(class_test_1)'
    xplot_1 = [xplot_1 distances_bw_1(class_test_1 == vis, axes(1))];
    yplot_1 = [yplot_1 distances_bw_1(class_test_1 == vis, axes(2))];
end

%second half
xplot_2 = []; yplot_2 = [];
for vis = unique(class_test_2)'
    xplot_2 = [xplot_2 distances_bw_2(class_test_2 == vis, axes(1))];
    yplot_2 = [yplot_2 distances_bw_2(class_test_2 == vis, axes(2))];
end


%plot all
for tw = 1:size(min([xplot_1 xplot_2]))
    for half = 1:2
       for vis = 1:4
           if half == 1
               plot(xplot_1(tw, vis), yplot_1(tw,vis), '.', 'markersize', mkrsize, 'linewidth', 1, 'color', colors(vis,:)) 
           else
               plot(xplot_2(tw, vis), yplot_2(tw,vis), '.', 'markersize', mkrsize, 'linewidth', 1, 'color', colors(vis,:)) 
           end
       end
    end
end

%
for vis = unique(ctx_context_bw)'
    for ic = unique(ctx_idx_bw)'
        %mean
        plot(mean(xplot_1(:, vis)), mean(yplot_1(:, vis)), '.', 'markersize', 3*mkrsize, 'linewidth', 5, 'color', colors(vis,:))
        plot(mean(xplot_1(:, vis)), mean(yplot_1(:, vis)), 'o', 'markersize', 3*mkrsize/pi, 'linewidth', 5, 'color', [0 0 0])
        
        plot(mean(xplot_2(:, vis)), mean(yplot_2(:, vis)), '.', 'markersize', 3*mkrsize, 'linewidth', 5, 'color', colors(vis,:))
        plot(mean(xplot_2(:, vis)), mean(yplot_2(:, vis)), 'o', 'markersize', 3*mkrsize/pi, 'linewidth', 5, 'color', [0 0 0])
    end
end
%}

plot([0 2], [0 2], 'k--', 'linewidth', 4)
hold on; plot([1 1], [0 2], 'k--', 'linewidth', 1)
hold on; plot([0 2], [1 1], 'k--', 'linewidth', 1)