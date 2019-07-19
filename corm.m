function cm = corm(c1, c2)
%function corm(cell_rates_first, cell_rates_second)
%
% makes correllation matrix of rates in c1 and c2

%preallocate
cm = nan(size(c1,2), size(c2,2));

%drop cells with unlikely firing rates
%{
min_rate = 1;
max_rate = 100;

min_idx = sum(c1<min_rate,2)>0 & sum(c2<min_rate,2)>0;
max_idx = sum(c1>max_rate,2)>0 & sum(c2>max_rate,2)>0;

c1 = c1(~min_idx & ~max_idx, :);
c2 = c2(~min_idx & ~max_idx, :);
%}

%normalize
%c1 = c1 - repmat(mean(c1, 2), 1, size(c1,2));%remove c1 means
%c1 = c1./repmat(std(c1, 0, 2), 1, size(c1,2)); %remove c1 stdevs
%c2 = c2 - repmat(mean(c2, 2), 1, size(c2,2));%remove c2 means
%c2 = c2./repmat(std(c2, 0, 2), 1, size(c2,2)); %remove c2 stdevs

%c = [c1 c2];
%c = c - repmat(mean(c1, 2), 1, size(c,2));%remove c1 means
%c = c./repmat(std(c, 0, 2), 1, size(c,2)); %remove c1 stdevs

%rate differences
%centerrate = mean(c(:));

%c1 = c(:,1:size(c1,2))
%c2 = c(:,(size(c1,2)+1):end)

%[R p] = fit_line(c1(:,1), c2(:,2)); title('c1 c2')
%[R p] = fit_line(c1(:,1), c2(:,3)); title('c1 c3')
%[R p] = fit_line(c1(:,2), c2(:,4)); title('c2 c4')

%correllate first and second half cell activity in each context
for context_first = 1:size(c1,2)
    
    for context_second = 1:size(c2,2)
        
        cm(context_first, context_second) = corr(c1(:,context_first), c2(:, context_second));
        
    end
    
end

mean([dist(cm(:,1)', cm(:,2)) dist(cm(:,3)', cm(:,4))])
mean([dist(cm(:,1)', cm(:,3)) dist(cm(:,2)', cm(:,4))])


%mirror
%cm = (cm + cm') ./2;

%plot
%
figure;
imagesc(cm)
colormap jet
colorbar
caxis([-1 1])
axis square
set(gca,'TickLength',[0, 0]);
%}

figure; 
%bar([mean([cm(1,2), cm(3,4)]) mean([cm(1,3), cm(1,4), cm(2,3), cm(2,4)])]) %redundant
bar([mean([cm(3,1), cm(4,2) cm(1,3), cm(2,4)]) mean([cm(2,1), cm(3,2), cm(4,3), cm(4,1) cm(1,2), cm(2,3), cm(3,4), cm(1,4)])]) %firsthalf-secondhalf BWBW
%bar([mean([cm(2,4) cm(4,2)]) mean([cm(1,3) cm(3,1)])]) %firsthalf-secondhalf AOAO

%{
recognize_arranges_ctxts = mean([cm(1,1) cm(3,3)]);
recognize_object_ctxts = mean([cm(2,2) cm(4,4)]);
discriminate_arrange_ctxts = mean([cm(1,3) cm(3,1)])
discriminate_object_ctxts = mean([cm(2,4) cm(4,2)])

all_recognize = mean([recognize_arranges_ctxts recognize_object_ctxts]);
double_discrim = mean([cm(1,2) cm(2,3) cm(3,4) cm(2,1) cm(3,2) cm(4,3) cm(1,4) cm(4,1)])

bar([discriminate_arrange_ctxts discriminate_object_ctxts double_discrim])
%bar([(recognize_arranges_ctxts - discriminate_arrange_ctxts) (recognize_object_ctxts-discriminate_object_ctxts) (all_recognize - double_discrim)])
%}
ylim([-1 1])

end