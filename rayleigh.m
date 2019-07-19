function [pval,R] = rayleigh(U)
%U is a vector length 360 or some multiple. Each item in the vector is a
%firing rate indicating the mean firing rate observed when the rat was
%facing that HD.




%find bin size assuming that the test is on 360 data
bin_size = 360/length(U);

%convert bins to radians
bins = deg2rad(1:bin_size:360);

%...and then to vectors
vects_x = cos(bins);
vects_y = sin(bins);

%fill matrix with a number of vects for each angle equal to the firing rate
vects = nan(sum(floor(U(:))),2);
iv_ct = 1;
for iv = 1: length(U)

    %floor(U(iv))
    
    vects(iv_ct:iv_ct+floor(U(iv))-1, :) = repmat([vects_x(iv) vects_y(iv)], floor(U(iv)), 1);
    iv_ct = iv_ct+floor(U(iv));
end

%resultant vector
r_vect = abs(mean(vects));

%resultant vector length
R = pdist([0 0; r_vect]);



%p value
pval = exp(sqrt(1+4*size(vects,1)+4*(size(vects,1)^2-R^2))-(1+2*size(vects,1)));
