function [peaks] = hd_tuning_new(rd)
%plots hd tuning of each cell

peaks = nan(length(rd),16); 


for icell = 1:161 

    %mean objHD tuning curves for each HD cell
    %{
    rd_obj_distrs = nanmean([cell2mat(rd{icell}(1:4,1));...
        cell2mat(rd{icell}(1:4,2));
        cell2mat(rd{icell}(1:4,3));
        cell2mat(rd{icell}(1:4,4))]);
    rd_obj_distrs = interp1(360/length(rd_obj_distrs):360/length(rd_obj_distrs):360, rd_obj_distrs, 1:360);
    rd_obj_distrs = smooth_around(rd_obj_distrs, 20);
    peaks(icell,1) = find(rd_obj_distrs==max(rd_obj_distrs));
    %}
close all

    %each object individually
    for ctx = 1:4
        for obj = 1:4
            rd_obj_distrs = cell2mat(rd{icell}(obj,ctx));
            rd_obj_distrs = interp1(360/length(rd_obj_distrs):360/length(rd_obj_distrs):360, rd_obj_distrs, 1:360);
            rd_obj_distrs = smooth_around(rd_obj_distrs, 20);
            peaks(icell, 4*(ctx-1)+obj) = find(rd_obj_distrs==max(rd_obj_distrs));
        end
    end
end
