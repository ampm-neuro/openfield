function [corms, same_diag_corr, withins_diag_corr, betweens_diag_corr] = all_pop_pixlecorr_plots(out_bw_first,out_bw_second, out_bw_all)

%standardize
out_bw_first = stdz_mtx(out_bw_first);
out_bw_second = stdz_mtx(out_bw_second);
out_bw_all = stdz_mtx(out_bw_all);


if nargin > 3
    shuf = varargin{1};
else
    shuf = 0;
end


%remove nan cells
%
all_mtx = sum(mean(cat(3, out_bw_first, out_bw_second, out_bw_all),3));
isntnan_idx = ~isnan(all_mtx); %for columns (cells)
find(isntnan_idx==0)

out_bw_first = out_bw_first(:, isntnan_idx, :);
out_bw_second = out_bw_second(:, isntnan_idx, :);
out_bw_all = out_bw_all(:, isntnan_idx, :);
%}

%BW

%withins
within_comps = [1 2; 3 4];
between_comps = [1 3; 1 4; 2 3; 2 4];

%preallocate correlation matrices
same_within = nan(size(out_bw_first,1), size(out_bw_first,1), size(out_bw_first,3));
withins1 = nan(size(out_bw_first,1), size(out_bw_first,1), size(within_comps,1));
betweens1 = nan(size(out_bw_first,1), size(out_bw_first,1), size(between_comps,1));
withins2 = nan(size(out_bw_first,1), size(out_bw_first,1), size(within_comps,1));
betweens2 = nan(size(out_bw_first,1), size(out_bw_first,1), size(between_comps,1));
withins = nan(size(out_bw_first,1), size(out_bw_first,1), size(within_comps,1));
betweens = nan(size(out_bw_first,1), size(out_bw_first,1), size(between_comps,1));


%compare sessions to themselves - how good is spatial decoding?
for ci = 1: size(out_bw_first,3)
    for si1 = 1:size(out_bw_first,1)
                
        %for si2 = 1:size(out_bw_first,1)
            same_within(si1, si1, ci) = corr(out_bw_first(si1,:,ci)', out_bw_second(si1,:,ci)');
        %end
    end
end


%plot 
comb_all_within = mean(same_within,3);
%{
figure; imagesc(comb_all_within); colormap jet; axis square
max_pos = nan(size(all_within,1),1);
for mpi = 1:size(comb_all_within,1)
    max_pos(mpi) = find(comb_all_within(mpi,:) == max(comb_all_within(mpi,:)), 1); %max col in row
end
hold on; plot(max_pos, 1:size(comb_all_within,1), 'k-', 'linewidth', 5)
%}
same_diag_corr = mean(comb_all_within(diag_mask(size(comb_all_within,2))));


%compare sessions to within and betweenosite type - can we decode across type?
%
%withins
for ci = 1: size(within_comps,1)
    for si = 1:size(out_bw_first,1)
        withins(si, si, ci) = corr(out_bw_all(si,:,within_comps(ci,1))', out_bw_all(si,:,within_comps(ci,2))'); 
    end
end
 comb_all_within = mean(withins,3);
 withins_diag_corr = mean(comb_all_within(diag_mask(size(comb_all_within,2))));
 
%betweenosites
for ci = 1: size(between_comps,1)
    for si = 1:size(out_bw_first,1)
        betweens(si, si, ci) = corr(out_bw_all(si,:,between_comps(ci,1))', out_bw_all(si,:,between_comps(ci,2))'); 
    end
end
 comb_all_between = mean(betweens,3);
 betweens_diag_corr = mean(comb_all_between(diag_mask(size(comb_all_between,2))));
%}



for ci = 1: size(within_comps,1)
    for si1 = 1:size(out_bw_first,1)
        %for si2 = 1:size(out_bw_first,1)
            withins1(si1, si1, ci) = corr(out_bw_first(si1,:,within_comps(ci,1))', out_bw_second(si1,:,within_comps(ci,2))');
            withins2(si1, si1, ci) = corr(out_bw_first(si1,:,within_comps(ci,2))', out_bw_second(si1,:,within_comps(ci,1))');
        %end
    end
end
 comb_all_within1 = mean(withins1,3);
 comb_all_within2 = mean(withins2,3);
 comb_all_within_within = mean(cat(3,comb_all_within1, comb_all_within2),3);
 %withins_diag_corr = mean(comb_all_within_within(diag_mask(size(comb_all_within_within,2))));


for ci = 1: size(between_comps,1)
    for si1 = 1:size(out_bw_first,1)
        %for si2 = 1:size(out_bw_first,1)
            betweens1(si1, si1, ci) = corr(out_bw_first(si1,:,between_comps(ci,1))', out_bw_second(si1,:,between_comps(ci,2))');
            betweens2(si1, si1, ci) = corr(out_bw_first(si1,:,between_comps(ci,2))', out_bw_second(si1,:,between_comps(ci,1))');
        %end
    end
end
 comb_all_within1 = mean(betweens1,3);
 comb_all_within2 = mean(betweens2,3);
 comb_all_within_between = mean(cat(3,comb_all_within1, comb_all_within2),3);
 %betweens_diag_corr = mean(comb_all_within_between(diag_mask(size(comb_all_within_between,2))));
 

%corm cell output
corms{1} = same_within;
corms{2} = withins1;
corms{3} = withins2;
corms{4} = betweens1;
corms{5} = betweens2;



%figures
%
%on-diag bar
%figure; bar([same_diag_corr, withins_diag_corr, betweens_diag_corr]);
figure; bar([withins_diag_corr, betweens_diag_corr]);

%same corm
%{
sames_mean = mean(same_within,3);
figure; imagesc(reshape(sames_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title same
%}

%withins corm
%{
withins_mean = mean(cat(3, withins1, withins2),3);
figure; imagesc(reshape(withins_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title withinhalves
%}
withins_mean = mean(withins,3);
figure; imagesc(reshape(withins_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title withinall

%betweens corm
%{
betweens_mean = mean(cat(3, betweens1, betweens2),3);
figure; imagesc(reshape(betweens_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title betweenhalves
%}
betweens_mean = mean(cat(3, betweens),3);
figure; imagesc(reshape(betweens_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title betweenall

%{
figure; imagesc(reshape(withins_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1)))-reshape(betweens_mean(diag_mask(size(out_bw_first,1))), sqrt(size(out_bw_first,1)), sqrt(size(out_bw_first,1))))
colormap jet; colorbar; caxis([-0.1 1]); axis square;
title subtract
%}


%subfunctions
    function mtxout = stdz_mtx(mtxin)
        mtxin = mtxin - mean(mtxin); 
        mtxout = mtxin./std(mtxin);
    end

end