function plot_peak_amp(peaks, rec_threshold, times)
fit_line(times, peaks);
hold on;
plot(xlim, [rec_threshold rec_threshold], '--')
plot(xlim, [min(peaks) min(peaks)], '--') %approx sorting threshold