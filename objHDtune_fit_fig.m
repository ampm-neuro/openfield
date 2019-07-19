function objHDtune_fit_fig(neuron, context, frame, rd, ...
    rs_all_obj, context_ia, all_counts_ia, all_pos_ia, all_hd_ia, all_obj_ia)

%ObjHD cells
objHD_idx = rs_all_obj(:,2)<0.05;

%mean objHD tuning curves for each HD cell
%
rd_obj = rd(objHD_idx);
rd_obj_distrs = (nanmean([cell2mat(rd_obj{neuron}(1:4, 5-4));...
    cell2mat(rd_obj{neuron}(1:4, 6-4));
    cell2mat(rd_obj{neuron}(1:4, 7-4));
    cell2mat(rd_obj{neuron}(1:4, 8-4))]));
rd_obj_distrs = interp1(360/length(rd_obj_distrs):360/length(rd_obj_distrs):360, rd_obj_distrs, 1:360);
rd_obj_distrs = smooth_around(rd_obj_distrs, 20);
figure; polarplot(norm_mtx(rd_obj_distrs'))


%relevant info at each time point
TW_FRs = all_counts_ia(context_ia==context, objHD_idx).*4; TW_FRs = TW_FRs(:,neuron);
TW_FRs = smooth(TW_FRs,10);
TW_posXY = all_pos_ia(context_ia==context, objHD_idx, :); TW_posXY = squeeze(TW_posXY(:,neuron,:));
TW_alloHD = all_hd_ia(context_ia==context, objHD_idx); TW_alloHD = TW_alloHD(:,neuron);
TW_objHD1234 = all_obj_ia(context_ia==context, 5:8, objHD_idx); TW_objHD1234 = TW_objHD1234(:,:,neuron);
TW_objdist1234 = all_obj_ia(context_ia==context, 1:4, objHD_idx); TW_objdist1234 = TW_objdist1234(:,:,neuron);

%stds
Std_FRs = std(TW_FRs);
    
%calculate expected rates for each object at each TW
angles_ihave = ((1:length(rd_obj_distrs)).*360/length(rd_obj_distrs)) - (360/length(rd_obj_distrs))/2;
rates_ihave = rd_obj_distrs;
expected_rates = nan(size(TW_objHD1234));
for itw = 1:size(TW_objHD1234,1)
    expected_rates(itw,:) = interp1(angles_ihave, rates_ihave, TW_objHD1234(itw,:));
end

for i = 1:size(expected_rates,2)
   expected_rates(:,i) = smooth(expected_rates(:,i),10);
end

%calculate abs(error) (distance between expected and observed rates)
error_FRs = abs(repmat(TW_FRs, 1, 4) - expected_rates);
error_FRs = error_FRs./Std_FRs;

%interpolate
current_time = .25:.25:length(TW_FRs)*.25;
desired_time = .25:.01:length(TW_FRs)*.25;


TW_FRs = interp1(current_time, TW_FRs, desired_time);
txy = TW_posXY; clearvars TW_posXY
TW_posXY(:,1) = interp1(current_time, txy(:,1), desired_time);
TW_posXY(:,2) = interp1(current_time, txy(:,2), desired_time);
TW_alloHD = circ_interp1(TW_alloHD', [0 360], desired_time', current_time');
error_FRs = interp1(current_time, error_FRs, desired_time);
er = expected_rates; clearvars expected_rates
expected_rates(:,1) = interp1(current_time, er(:,1), desired_time);
expected_rates(:,2) = interp1(current_time, er(:,2), desired_time);
expected_rates(:,3) = interp1(current_time, er(:,3), desired_time);
expected_rates(:,4) = interp1(current_time, er(:,4), desired_time);


%norm rates
FR_norm = norm_mtx(TW_FRs(1:end-100));
ER_norm = norm_mtx(expected_rates);
ER_sumnorm = mean(ER_norm,2);

figure;
rng = [-(range(frame)-range(frame)/2) (range(frame)-range(frame)/2)].*.01;
colors = get(gca,'ColorOrder');
subplot(1,2, 1)
object_grid(context);
set(gca,'TickLength',[0, 0]); box off
hold on
subplot(1,2, 2)
set(gca,'TickLength',[0, 0]); box off; ylim([0 1]); 

for TW = frame
    subplot(1,2, 1)
hold on
    angle = TW_alloHD(TW);
    s = TW_posXY(TW,1);
    t = TW_posXY(TW,2);
    u = sin(angle*(pi/180));
    v = cos(angle*(pi/180));
    line_length = .09;
    arrow_size = 0;
    
    %plot pos and hd
    h1 = quiver(s, t, u, v, line_length, 'k', 'linewidth', 3);
    quiver_arrow_size(h1, arrow_size)
    text(s, t, num2str(find(frame==TW)-1))
    
    %bars
    if TW == median(frame)
    sp2 = subplot(1,2, 2);
    %ylim([0 3]); xlim([0 5])
    
    local_error = error_FRs(TW,:);
    hold on
    %for i = 1:length(local_error)
    %    bar(i, local_error(i), 'facecolor', colors(i,:));
    %end
    if TW>50
        plot_rng = rng(1)/.01:1:rng(2)/.01;
        plot(plot_rng, FR_norm((TW+rng(1)/.01):1:(TW+rng(2)/.01)), 'k-')
        plot(plot_rng, ER_sumnorm((TW+rng(1)/.01):1:(TW+rng(2)/.01)), '-', 'color', [.5 .5 .5])
        %plot([0 0], ylim, 'k-', 'linewidth', 3)
        for i = 1:length(local_error)
            plot(plot_rng, ER_norm((TW+rng(1)/.01):1:(TW+rng(2)/.01), i), '-', 'color', colors(i,:))
        end
    end
    hold off
    end

end

end