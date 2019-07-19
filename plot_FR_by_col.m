function plot_FR_by_col(eptrials, cluster, col, bins, tw_size, slide_size, zscore_yn)
%plot FR_by_col for one cell, overlaying each context visit
%
%   tw and slide should be in ms. 500 is 500ms or .5s
%
%   zscore_yn determines whether rates are zscored. 1 for z, 0 for no.
%
%   eptrials(:,7) contains the rat's velocity.
%   eptrials(:,8) contains the rat's acceleration
%   eptrials(:,9) contains the rat's head direction.

%calculate binned rates
rates_all = [];
col_all = [];
ctxt_idx = [];
means_all = nan(4, 1);
stds_all = nan(4, 1);
for ctxt = 1:4

    %bin FRs by col
    cur_eptrials = eptrials(eptrials(:,5)==ctxt, :);
    [rates_out, col_out] = FR_by_col(cur_eptrials, cluster, col, bins, tw_size, slide_size);

    %save means and stds for later
    means_all(ctxt, :) = nanmean(rates_out, 1);
    stds_all(ctxt, :) = nanstd(rates_out, 1);

    %zscore each context independently
    if zscore_yn == 1
        rates_out = zscore_mtx(rates_out);
    elseif zscore_yn == 0
    else
        error('zscore_yn must be 1 to zscore or 0 to not.')
    end
    
    %combine normed contexts
    rates_all = [rates_all; rates_out];
    col_all = [col_all; col_out];
    ctxt_idx = [ctxt_idx; repmat(ctxt, size(col_out))];
    
end

%absolute value
%col_all = abs(col_all);

%plot
if col == 7
    degree = 2;
elseif col == 8
    degree = 4;
end
fit_line(col_all, rates_all, degree);

%add context colored dots
hold on;
set(gca,'TickLength',[0, 0]);
for ctxt = 1:4
    plot(col_all(ctxt_idx==ctxt), rates_all(ctxt_idx==ctxt), '.', 'markersize', 15)
end
set(gca,'TickLength',[0, 0]);