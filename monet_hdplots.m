function monet_hdplots(all_pos_ia, all_hd_ia, context_ia, context, neuron)
%plots path and overlays colored dots corresponding to HD
figure; hold on; 
plot(all_pos_ia(context_ia==context,neuron,1), all_pos_ia(context_ia==context,neuron,2)...
    ,'-', 'color', [.9 .9 .9])

hd_rng = 0:10:360; 
for i = 2:length(hd_rng)
    idx = context_ia==context & all_hd_ia(:,neuron)>hd_rng(i-1)...
        & all_hd_ia(:,neuron)<hd_rng(i);
    plot(all_pos_ia(idx,neuron,1), all_pos_ia(idx,neuron,2), '.', 'color',...
        [[1 1].*(hd_rng(i)/360) .5], 'markersize', 30)
end
axis square off; %set(gca,'TickLength',[0, 0]); box off