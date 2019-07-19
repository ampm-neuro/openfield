function [rotation, rot_cor, orig_cor] = rot_fit(rd)
%takes 2 rate distriubtions and finds the clockwise angular rotation (of 
%distribution 2) that maximizes the correllation between them

basic_idx = 1:size(rd,2);
CW_rot_deg = 0:359;
corrs = nan(length(CW_rot_deg),1);

%rotate and correlate
for ir = CW_rot_deg
    
   %rotate
   rd2 = rd(2,circshift(basic_idx,ir));
   
   %correlate
   corrs(ir+1) = corr(rd(1,:)', rd2');
    
end

rotation = CW_rot_deg(corrs==max(corrs));
rot_cor = corrs(corrs==max(corrs));
orig_cor = corrs(1);

%{
figure; hold on
plot(rd(1,:))
plot(rd(2,:))
plot(rd(2,circshift(basic_idx,rotation)))
%}
end