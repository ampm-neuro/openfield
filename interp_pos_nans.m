function eptrials = interp_pos_nans(eptrials)
%interpolate_at_nans finds missing position information in eptrials
%and interpolates to fill
%

    %index for missing position values
    real_x = eptrials(~isnan(eptrials(:,2)),2);
    real_y = eptrials(~isnan(eptrials(:,3)),3);
    all_time = eptrials(:,1);

    %only including unique time points among the missing rows
    [real_time_x, idx_x] = unique(eptrials(~isnan(eptrials(:,2)),1));
    [real_time_y, idx_y] = unique(eptrials(~isnan(eptrials(:,3)),1));

    %index real_x and real_y to match real_time unique elements
    eptrials(:,2) = interp1(real_time_x, real_x(idx_x), all_time);
    eptrials(:,3) = interp1(real_time_y, real_y(idx_y), all_time);  

end