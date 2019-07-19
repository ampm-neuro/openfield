function [spatial_corrs, EFRs, objangs] = objHDtune_space(neuron, rd_in, rs_all_obj_in, rate_matrices_ia_in, hd_proportions_ia_in)
%calculate objHD tuning curve for ObjHD neuron (input)
%use to predict spatial maps context (input)

%load('ampm_HDstuff_2_18_18.mat')
%load('250ms_withallobj.mat')

spatial_corrs = nan(4,1);

figure;
for ctx_count = 1:4
    context = ctx_count+4;
    
    %reset vars
    rd = rd_in;
    rs_all_obj = rs_all_obj_in; 
    rate_matrices_ia = rate_matrices_ia_in; 
    hd_proportions_ia = hd_proportions_ia_in;
    
objHD_idx = logical(ones(size(rs_all_obj(:,2)<0.05)));
rate_matrices_ia = rate_matrices_ia(:,:,objHD_idx);

%mean objHD tuning curves for each HD cell
%
rd_obj = rd(objHD_idx);
rd_obj_distrs = ([cell2mat(rd_obj{neuron}(1:4, 5-4));...
    cell2mat(rd_obj{neuron}(1:4, 6-4));
    cell2mat(rd_obj{neuron}(1:4, 7-4));
    cell2mat(rd_obj{neuron}(1:4, 8-4))]);
rd_obj_distrs = norm_mtx(rd_obj_distrs')';
rd_obj_distrs = nanmean(rd_obj_distrs(setdiff(5:8,context),:));
rd_obj_distrs = norm_mtx(rd_obj_distrs')';
rd_obj_distrs = interp1(360/length(rd_obj_distrs):360/length(rd_obj_distrs):360, rd_obj_distrs, 1:360);
rd_obj_distrs = smooth_around(rd_obj_distrs, 20);

%HD proportions (becomes bins*bins*360)
hd_proportions_ia = squeeze(hd_proportions_ia(context-4,:,:,:,objHD_idx));
hd_proportions_ia = hd_proportions_ia(:,:,:,neuron);

%hard-coded average object locations for all sessions
object_locations(1,:,1) = [mean([0.5 0.6875]) mean([0.6875 1-0.125])] + [-0.005 0];
object_locations(2,:,1) = [mean([0.3125 0.5]) mean([0.5 0.6875])] + [-0.005 0.005];
object_locations(3,:,1) = [mean([0.6875 1-0.125]) mean([0.3125 0.5])] + [0 0.005];
object_locations(4,:,1) = [mean([0.3125 0.5]) mean([0.125 0.3125])] + [-0.005 0.005];
object_locations(1,:,2) = [mean([0.5 0.6875]) mean([0.5 0.6875])] + [-0.0025 0.0025];
object_locations(2,:,2) = [mean([0.6875 1-0.125]) mean([0.5 0.6875])] + [0 0.005];
object_locations(3,:,2) = [mean([0.125 0.3125]) mean([0.3125 0.5])] + [-0.0075 0.005];
object_locations(4,:,2) = [mean([0.3125 0.5]) mean([0.3125 0.5])] + [0 0.005];
object_locations(1,:,3) = [mean([0.125 0.3125]) mean([0.6875 1-0.125])] + [-0.0025 0.005];
object_locations(2,:,3) = [mean([0.6875 1-0.125]) mean([0.6875 1-0.125])] + [-0.0075 0.0075];
object_locations(3,:,3) = [mean([0.125 0.3125]) mean([0.125 0.3125])] + [-0.0075 0.0025];
object_locations(4,:,3) = [mean([0.6875 1-0.125]) mean([0.125 0.3125])] + [-0.005 -0.005];
object_locations(1,:,4) = [mean([0.125 0.3125]) mean([0.6875 1-0.125])] + [0 0.005];
object_locations(2,:,4) = [mean([0.6875 1-0.125]) mean([0.6875 1-0.125])] + [-0.0075 0.0075];
object_locations(3,:,4) = [mean([0.125 0.3125]) mean([0.125 0.3125])] + [-0.005 0];
object_locations(4,:,4) = [mean([0.6875 1-0.125]) mean([0.125 0.3125])] + [-0.005 -0.005];

%position range
pos_rng = linspace(0,1,20);

%object angles
objangs = nan(length(pos_rng), length(pos_rng), 4, 360);

%iterate through positions
rat_ang = 1;
for rat_pos_x = pos_rng
    for rat_pos_y = pos_rng

        rat_pos = [rat_pos_x, rat_pos_y];

        %all object positions in this context
        object_pos = object_locations(:,:,context-4)';   
        %
        if sum(pdist([rat_pos; object_locations(:,:,context-4)])<0)>0
            objangs(pos_rng==rat_pos_y, pos_rng==rat_pos_x, :, rat_ang) = nan;
            continue
        end
        %}

        %position difference to center of objects
        pos_difference = repmat(reshape(object_pos,1,2,4), ...
            size(rat_pos,1), 1, 1) - repmat(rat_pos, 1, 1, 4); %x and y offset between rat and object

        %angle between rat facing north (0deg) and object
        objangs(pos_rng==rat_pos_y, pos_rng==rat_pos_x, :, rat_ang) = ...
            atan2d(squeeze(pos_difference(:, 1, :)), squeeze(pos_difference(:, 2, :))); 
    end
end
for rat_ang = 2:360
    objangs(:,:,:,rat_ang) = objangs(:,:,:,rat_ang-1) - 1;
end
objangs(objangs<0) = rem(objangs(objangs<0),360) + 360; 
objangs(objangs<1) = round(objangs(objangs<1));
objangs(objangs>360.5) = 1; 
objangs(objangs>360 & objangs<360.5) = 360;
objangs(objangs==0) = 360;

%calculate expected firing rates
%4d matrix (ybin,xbin, obj_page, rat_angle) 
% each page is a matrix of expected rates at all positions for 1 object
% four pages (1 per object)
% 4th dim is 1:360 deg of the rat's current HD
EFRs = interp1(1:360, rd_obj_distrs, objangs);
EFRs = mean(EFRs,3); 
EFRs = squeeze(EFRs); %bins,bins,360


%weight by the distribution of HDs in each spatial bin
EFRs_mean = nan(size(EFRs,1), size(EFRs,1));
for iy = 1:size(EFRs,1)
    for ix = 1:size(EFRs,1)
        observed_hd_idx = ~isnan(hd_proportions_ia(iy,ix,:));
        hd_proportions = hd_proportions_ia(iy,ix,observed_hd_idx)...
            ./nansum(hd_proportions_ia(iy,ix,:));
        EFRs_at_ohds = EFRs(iy,ix, observed_hd_idx);
        EFR_weightedmean = sum(EFRs_at_ohds.*hd_proportions);
        EFRs_mean(iy,ix) = EFR_weightedmean;
    end
end
EFRs = smooth2a(EFRs_mean,1);
EFRs = inpaint_nans(EFRs);
EFRs = reshape( norm_mtx(EFRs(:)), size(EFRs,1), size(EFRs, 2));

%plot side by side
%figure; hold on
%true heatmap
subplot(4,2,2*(ctx_count-1)+1)
imagesc(EFRs)
object_grid(context,size(EFRs,1));
caxis([0 1])
axis square; set(gca,'TickLength',[0, 0]); axis off
title('Predicted')

%HD estimate
subplot(4,2,2*(ctx_count-1)+2)
sidelgth = sqrt(length(rate_matrices_ia(context-4, :, neuron)));
rm_obs = reshape(rate_matrices_ia(context-4, :, neuron), sidelgth, sidelgth);
rm_obs = smooth2a(rm_obs,1);
%rm_obs = inpaint_nans(rm_obs);
rm_obs = flipud(reshape( norm_mtx(rm_obs(:)), size(rm_obs,1), size(rm_obs, 2)));
imagesc(rm_obs)
%object_grid(context,size(rm_obs,1));
caxis([0 1])
axis square; set(gca,'TickLength',[0, 0]); axis off
title('Observed')

EFRs_corr_in = EFRs(:);
rm_obs_corr_in = rm_obs(:);
nnan_idx = ~isnan(EFRs_corr_in) & ~isnan(rm_obs_corr_in);
spatial_corrs(ctx_count) = corr(EFRs_corr_in(nnan_idx), rm_obs(nnan_idx));

end
set(gcf, 'Position', [600, 200, 400, 800])




end






