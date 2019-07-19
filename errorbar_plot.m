function errorbar_plot( cell_e )
%plots errorbar

%preallocate
all_mtx = cell2mat(cell_e(:));
all_mtx_grps = [];

means = nan(length(cell_e),1);
stds = nan(length(cell_e),1);
sqrt_l = nan(length(cell_e),1);


%figure; 
hold on;

for imtx = 1: length(cell_e)
  
    for i2 = 1:length(cell_e{imtx})
        plot(imtx+(rand(1)/4-rand(1)/4), cell_e{imtx}(i2), 'ko')
    end
    
    all_mtx_grps = [all_mtx_grps; repmat(imtx, size(cell_e{imtx}(:)))];
    means(imtx) = nanmean(cell_e{imtx});
    stds(imtx) = nanstd(cell_e{imtx});
    sqrt_l(imtx) = sqrt(length(~isnan(cell_e{imtx})));

    
end

std_es = stds./sqrt_l;


set(gca,'TickLength',[0, 0]); box off
xlim([.5 length(cell_e)+.5])
xticks(1:4)


    
    %errorbar(means, std_es, 'k')

%[P,ANOVATAB] = anova1(all_mtx,all_mtx_grps, 'off')
%[p] = anova_rm(reshape(all_mtx, length(all_mtx)/length(cell_e), length(cell_e)));


end

