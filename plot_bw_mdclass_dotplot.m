mkrsize = 30;
colors = ([...
0.15 0.15 0.15;
0.3 0.3 0.3;
0.75 0.75 0.75;
0.9 0.9 0.9;
]);

%FIRST TWO AXES
figure; hold on; axis([.6 1.4 .6 1.4]); axis square; set(gca,'TickLength',[0, 0]); box off
xlabel('distance to black 1')
ylabel('distance to white 1')
%each time window
axes = [1 2];
for itw = 1:length(ctx_idx_bw)/4
    %each context
    for ic = 1:4
        %hold
        xplot = distances_bw(ctx_idx_bw == ic, axes(1));
        yplot = distances_bw(ctx_idx_bw == ic, axes(2));
        %plot        
        plot(xplot(itw), yplot(itw), '.', 'markersize', mkrsize, 'linewidth', 1, 'color', colors(ic,:))
        %plot(xplot(itw), yplot(itw), 'o', 'markersize', mkrsize/pi, 'linewidth', 1, 'color', [0 0 0])
    end
end

for ic = 1:4
    %mean
    plot(mean(distances_bw(ctx_idx_bw==ic, axes(1))), mean(distances_bw(ctx_idx_bw==ic, axes(2))), '.', 'markersize', 3*mkrsize, 'linewidth', 5, 'color', colors(ic,:))
    plot(mean(distances_bw(ctx_idx_bw==ic, axes(1))), mean(distances_bw(ctx_idx_bw==ic, axes(2))), 'o', 'markersize', 3*mkrsize/pi, 'linewidth', 5, 'color', [0 0 0])
end

plot([0 2], [0 2], 'k--', 'linewidth', 4)
hold on; plot([1 1], [0 2], 'k--', 'linewidth', 1)
hold on; plot([0 2], [1 1], 'k--', 'linewidth', 1)


%SECOND TWO AXES
figure; hold on; axis([.6 1.4 .6 1.4]); axis square; set(gca,'TickLength',[0, 0]); box off
xlabel('distance to black 2')
ylabel('distance to white 2')
axes = [3 4];
for itw = 1:length(ctx_idx_bw)/4
    %each timewindow
    for ic = 1:4
        %hold
        xplot = distances_bw(ctx_idx_bw == ic, axes(1));
        yplot = distances_bw(ctx_idx_bw == ic, axes(2));
        %plot        
        plot(xplot(itw), yplot(itw), '.', 'markersize', mkrsize, 'linewidth', 1, 'color', colors(ic,:))
        %plot(xplot(itw), yplot(itw), 'o', 'markersize', mkrsize/pi, 'linewidth', 1, 'color', [0 0 0])
    end
end

for ic = 1:4
    %mean
    plot(mean(distances_bw(ctx_idx_bw==ic, axes(1))), mean(distances_bw(ctx_idx_bw==ic, axes(2))), '.', 'markersize', 3*mkrsize, 'linewidth', 5, 'color', colors(ic,:))
    plot(mean(distances_bw(ctx_idx_bw==ic, axes(1))), mean(distances_bw(ctx_idx_bw==ic, axes(2))), 'o', 'markersize', 3*mkrsize/pi, 'linewidth', 5, 'color', [0 0 0])
end

plot([0 2], [0 2], 'k--', 'linewidth', 4)
hold on; plot([1 1], [0 2], 'k--', 'linewidth', 1)
hold on; plot([0 2], [1 1], 'k--', 'linewidth', 1)