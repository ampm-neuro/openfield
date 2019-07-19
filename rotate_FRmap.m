function [peak_r_rots, rvals] = rotate_FRmap(rm)
%average within contexts, carve a circle out of the rate matrices, rotate
%clockwise and take correlation each time. report rotation that produces
%highest spatial correlation

%smooth factor
sf = 3;

%circle mask
side_length = sqrt(size(rm,2));
radius = floor(side_length/2);
cm = createCirclesMask([side_length,side_length],...
    [side_length/2 + .5, side_length/2 + .5], radius-1);

%preallocate 
peak_r_rots = nan(size(rm,3),1); %peak rotation angles
within_spacecor = nan(size(rm,3),1);
between_spacecor = nan(size(rm,3),1);
between_rot_spacecor = nan(size(rm,3),1);

%iterate through cells
for mtx_num = 1:size(rm,3)
    
    %space corrs
    %
    b1_rm = reshape(rm(1, :, mtx_num), side_length, side_length); 
        b1_rm(~cm) = nan; b1_rm = smooth2a(b1_rm,sf);
    b2_rm = reshape(rm(2, :, mtx_num), side_length, side_length);
        b2_rm(~cm) = nan; b2_rm = smooth2a(b2_rm,sf);
    w1_rm = reshape(rm(3, :, mtx_num), side_length, side_length);
        w1_rm(~cm) = nan; w1_rm = smooth2a(w1_rm,sf);
    w2_rm = reshape(rm(4, :, mtx_num), side_length, side_length);
        w2_rm(~cm) = nan; w2_rm = smooth2a(w2_rm,sf);
    %within
    within1 = corrcoef(b1_rm(:),b2_rm(:),'rows','complete'); within1 = within1(2,1);
    within2 = corrcoef(w1_rm(:),w2_rm(:),'rows','complete'); within2 = within2(2,1);
    within_spacecor(mtx_num) = mean([within1 within2]);
    %between
    between1 = corrcoef(b1_rm(:),w1_rm(:),'rows','complete'); between1 = between1(2,1);
    between2 = corrcoef(b1_rm(:),w2_rm(:),'rows','complete'); between2 = between2(2,1);
    between3 = corrcoef(b2_rm(:),w1_rm(:),'rows','complete'); between3 = between3(2,1);
    between4 = corrcoef(b2_rm(:),w2_rm(:),'rows','complete'); between4 = between4(2,1);
    between_spacecor(mtx_num) = mean([between1 between2 between3 between4]);
    
    %average within context
    b_rm = nanmean(rm(1:2, :, mtx_num)); 
        b_rm = reshape(b_rm, side_length, side_length);
    w_rm = nanmean(rm(3:4, :, mtx_num));
        w_rm = reshape(w_rm, side_length, side_length);

    %smooth
    b_rm = smooth2a(b_rm,sf);
    w_rm = smooth2a(w_rm,sf);
    
    
    
    %preallocate r_vals
    rs = nan(360,1);

    %rotate and load r val
    for rot_angle = 0:359
        
        %rotate matrix
        rot_w_rm = rotate_about_center(w_rm, rot_angle);

        %carve out center circle
        b_rm(~cm) = nan;
        rot_w_rm(~cm) = nan;

        %correlate
        b_vect = b_rm(:);
        w_vect = rot_w_rm(:);
        nnan_idx = ~isnan(b_vect) & ~isnan(w_vect);
        rs(rot_angle+1) = corr(b_vect(nnan_idx), w_vect(nnan_idx));
    end
    
    %load rotation
    peak_rotangle = find(rs==max(rs)) -1;
    
    if numel(peak_rotangle)>1
        circdists = rad2deg(circ_dist(deg2rad([1 2 3]), deg2rad(359))); 
        circdists(circdists<360) = circdists(circdists<360)+360;
        peak_rotangle = peak_rotangle(circdists==min(circdists));
    end
    
    peak_r_rots(mtx_num) = peak_rotangle;
    
    
    
    
    %calculate betweens after rotating white contexts
    w1_rm = rotate_about_center(reshape(rm(3, :, mtx_num), side_length, side_length), peak_rotangle); 
        w1_rm(~cm) = nan; w1_rm = smooth2a(w1_rm,sf);
    w2_rm = rotate_about_center(reshape(rm(4, :, mtx_num), side_length, side_length), peak_rotangle); 
        w2_rm(~cm) = nan; w2_rm = smooth2a(w2_rm,sf);
    %between
    between1 = corrcoef(b1_rm(:),w1_rm(:),'rows','complete'); between1 = between1(2,1);
    between2 = corrcoef(b1_rm(:),w2_rm(:),'rows','complete'); between2 = between2(2,1);
    between3 = corrcoef(b2_rm(:),w1_rm(:),'rows','complete'); between3 = between3(2,1);
    between4 = corrcoef(b2_rm(:),w2_rm(:),'rows','complete'); between4 = between4(2,1);
    between_rot_spacecor(mtx_num) = mean([between1 between2 between3 between4]);
        
    
    if between_spacecor(mtx_num) < .2 && between_rot_spacecor(mtx_num)>.3
        %plot best fit
    	rot_fig(b_rm, w_rm, peak_rotangle, cm, sf)
    end
    
end

%combine space correllation r vals
rvals = [within_spacecor between_spacecor between_rot_spacecor];
end

function rot_fig(b_rm, w_rm, rot_angle, cm, sf)
    %plot
    rot_w_rm = rotate_about_center(w_rm, rot_angle);
    
    b_rm_hold = inpaint_nans(b_rm);
    b_rm_hold(~cm) = nan;
    w_rm_hold = inpaint_nans(w_rm);
    w_rm_hold(~cm) = nan;
    rot_w_rm_hold = inpaint_nans(rot_w_rm);
    rot_w_rm_hold(~cm) = nan;
    
    b_img = ones(size(b_rm_hold));
    b_img(isnan(b_rm_hold)) = 0;
    w_img = ones(size(w_rm_hold));
    w_img(isnan(w_rm_hold))=0;
    rot_w_img = ones(size(rot_w_rm_hold));
    rot_w_img(isnan(rot_w_rm_hold))=0;
    
    
    figure; 
    
    subplot(1,3,2)
    b_rm_hold = smooth2a(b_rm_hold,sf);
    imagesc(b_rm_hold, 'AlphaData', b_img)
    axis square off
    title('black')
    
    subplot(1,3,1)
    w_rm_hold = smooth2a(w_rm_hold,sf);
    imagesc(w_rm_hold, 'AlphaData', w_img)
    axis square off
    title('white')
    
    subplot(1,3,3)
    rot_w_rm_hold = smooth2a(rot_w_rm_hold,sf);
    imagesc(rot_w_rm_hold, 'AlphaData', rot_w_img)
    axis square off
    title(['white rot -' num2str(rot_angle) ' deg'])

end


function rm_out = rotate_about_center(rm_in, rot_angle)
%rotate the matrix about the center

% i j coordinates of matrix of rm_in size
ijs = sortrows([pairwise(1:size(rm_in,1)); ...
    [(1:size(rm_in,1))' (1:size(rm_in,1))']]);

% i j coordinates rotates by angle
ijs_rot = rotate_pts(rot_angle, ijs, [size(rm_in,1)/2+.5 size(rm_in,2)/2+.5]);

% round
ijs_rot = round(ijs_rot);

% index rotated i j coordinates that fall within bounds of rm_in
incl_ijs = ijs_rot(:,1) >= 1 & ijs_rot(:,1) <= size(rm_in,1)...
    & ijs_rot(:,2) >= 1 & ijs_rot(:,2) <= size(rm_in,2);

%preallocate new matrix
rm_out = nan(size(rm_in));

%iterate to load rm_out
ijs_ct = 1:size(ijs,1);
ijs_ct = ijs_ct(incl_ijs);
for i = ijs_ct
    rm_out(ijs_rot(i,1), ijs_rot(i,2)) = rm_in(ijs(i,1), ijs(i,2));
end



end