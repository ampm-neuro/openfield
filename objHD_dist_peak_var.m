[all_rd_heat] = ALL_combined_distHDobj;
all_var = nan(size(all_rd_heat,3),1);
for i = 1:size(all_rd_heat,3)
    mtx = all_rd_heat(:,:,i);
    %mtx = smooth2a(mtx,3)';
    %figure; hold on; imagesc(mtx);
    var_hold = nan(30, 1); 
    for i2 = 1:length(var_hold)
        find_hold = find(mtx(i2,:)==max(mtx(i2,:)), 1);
        var_hold(i2) = find(mtx(i2,:)==max(mtx(i2,:)), 1);
        %plot(find_hold, i2, 'ko', 'markersize', 20);
    end
    all_var(i) = circ_std(deg2rad(var_hold.*(360/size(mtx,2))));
    %title([num2str(all_var(i)) ' ' num2str(rs_all_obj(i,2))])
end
c{1} = rad2deg(all_var(rs_all_obj(:,2)>0.05)); 
c{2} = rad2deg(all_var(rs_all_obj(:,2)<=0.05)); 
c{3} = rad2deg(all_var(rs_all_allo(:,2)<=0.05));
errorbar_barplot(c([2 3 1])); errorbar_plot(c([2 3 1]))