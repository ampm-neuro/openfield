function arcFRsplot(FRs)

    figure; 
    fit_line_ratetime(FRs, 2);
    set(gcf,'pos',[0 600 700 500])
    set(gca,'TickLength',[0, 0]); box off

    %title([num2str(session_file),'-tt-' ,num2str(i),'-cell-' ,num2str(each_cell)])
    xlabel('Time from Start (min)')
    ylabel('Firing Rate (Hz)')

    %events
    hold on
    %for ev = 1:length(TimeStampsEV)
    %    plot([TimeStampsEV(ev) TimeStampsEV(ev)], ylim, 'r-')
    %end
    hold off
    xticklabels((xticks)/60)

end