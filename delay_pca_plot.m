function [all_mtx_z, EXPLAINED] = delay_pca_plot(all_time_windows, variables, pc_axes, bins)
%plots output from ALL_spikects_delay
%
% green is visited-left-reward
% blue is visited-right-reward
%
% over the delay, colors go from dark to light

%~input variables
num_intermediate_means = bins; %number of along-the-way means to plot
num_seconds = 0; %seconds to take off end(s) of delay
small_dot_size = 10; %scatterplot dot sizes for timewindow dots...
medium_dot_size = 50; %intermediate means and...
large_dot_size = 100; %overall means

%time_saving calculations
num_axes = length(pc_axes);
num_vars = length(variables);
axes_idx = 1:length(pc_axes);

%timewindow count
bins = length(all_time_windows(:,1,1));
bookend_bins = floor((num_seconds/30)*bins);

%cut off bookend_bins 
%all_time_windows = all_time_windows((bookend_bins+1):(end-bookend_bins), :, :); %off begining AND end
all_time_windows = all_time_windows(1:(end-bookend_bins), :, :); %off end only

%combine desired axes
all_mtx = [];
for ivar = 1:length(variables)
    all_mtx = [all_mtx; all_time_windows(:,:,variables(ivar))];
end

%zscore
all_mtx_z = all_mtx - repmat(mean(all_mtx), size(all_mtx,1),1);
all_mtx_z = all_mtx_z./std(all_mtx_z);

%reshape????
all_mtx_z = permute(reshape(all_mtx_z', size(all_time_windows,2), size(all_time_windows,1), length(axes)),[2 1 3]);

%pca
[all_mtx_pca, ~, ~, ~, EXPLAINED] = pca(all_mtx_z');

%preallocate imean holders
icell = cell(num_intermediate_means, num_vars);
imean = nan(num_intermediate_means, num_axes, num_vars);
icolor_hold = nan(num_intermediate_means, 3, num_vars);

%find rows corresponding to pre_zscored variables
num_wdw_all = size(all_mtx_pca, 1);
var_mtx = nan(num_wdw_all/num_vars, num_vars, num_axes);
idx_mtx = reshape(1:num_wdw_all, num_wdw_all/num_vars, num_vars);
for nvi = 1:num_vars
    var_mtx(:, nvi, :) = all_mtx_pca(idx_mtx(:,nvi), pc_axes);
end


%figure; imagesc(squeeze(var_mtx(:, 1, :)))
%figure; imagesc(squeeze(var_mtx(:, 2, :)))

%plot
figure; hold on;

%color modifier
base = [36 39 127;... %left correct
        17 127 35;... %right correct
        36 39 127;... %left error
        17 127 35]./300; %right error
    
icolor = .55/size(var_mtx,1);

%

%plots by number of dimensions

%iterate through vars
for i_v = 1:size(var_mtx,2)

    %plot time_windows
    %iterate over time (darker colors at begining)
    im_it = 1; im_ctr = 0; %counter for intermediate mean plots
    im_cc = base; %imean plot colors
    for i_t = 1:size(var_mtx,1)

    color_code = base(i_v,:) + repmat(icolor,1,3).*i_t;
        
        %use plot if 2d
        if num_axes==2


            %plot TIME WINDOWS
            %
            plot(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)),...
                '.', 'Color', color_code, 'MarkerSize', small_dot_size)

                %encircle dots
                if ismember(variables(i_v), [3 4]) %red for errors
                    plot(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)),...
                        'ro', 'MarkerSize', small_dot_size/3.5, 'LineWidth', .25)
                else %black for corrects
                    plot(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)),...
                        'ko', 'MarkerSize', small_dot_size/3.5, 'LineWidth', .25)
                end
            %}
 
            %plot INTERMEDIATE MEANS
            im_rng_hi = floor(size(var_mtx,1)/num_intermediate_means);            
            if ismember(i_t, floor(im_rng_hi:im_rng_hi:size(var_mtx,1)))
                im_ctr = im_ctr+1;
                
                %timewindow bins in intermediate ranges
                im_range = im_it:i_t;
                
                %load imean holders
                imean(im_ctr, :, i_v) = [mean(var_mtx(im_range, i_v, axes_idx(1))) mean(var_mtx(im_range, i_v, axes_idx(2)))];
                icolor_hold(im_ctr, :, i_v) = mean([color_code; im_cc(i_v,:)]);
                
                %plot imean dot
                plot(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))),...
                    '.', 'Color', icolor_hold(im_ctr, :, i_v), 'MarkerSize', medium_dot_size)
                
                    %encircle imean dot
                    if ismember(variables(i_v), [3 4]) %red for errors
                        plot(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))),...
                            'ro', 'MarkerSize', medium_dot_size/3.5, 'LineWidth', .75)
                    else %black for corrects
                        plot(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))),...
                            'ko', 'MarkerSize', medium_dot_size/3.5, 'LineWidth', .75)
                    end
                %reset counter
                im_it = i_t+1;
                im_cc(i_v,:) = color_code;
            end

        %use plot3 for 3d plots    
        elseif num_axes==3

            %plot TIME WINDOWS
            %
            plot3(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)), var_mtx(i_t, i_v, axes_idx(3)),...
                '.', 'Color', color_code, 'MarkerSize', small_dot_size)

                %encircle dots
                if ismember(variables(i_v), [3 4]) %red for errors
                    plot3(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)), var_mtx(i_t, i_v, axes_idx(3)),...
                        'ro', 'MarkerSize', small_dot_size/3.5, 'LineWidth', .25)
                else %black for corrects
                    plot3(var_mtx(i_t, i_v, axes_idx(1)), var_mtx(i_t, i_v, axes_idx(2)), var_mtx(i_t, i_v, axes_idx(3)),...
                        'ko', 'MarkerSize', small_dot_size/3.5, 'LineWidth', .25)
                end
            %}
            
            
            %plot INTERMEDIATE MEANS
            im_rng_hi = floor(size(var_mtx,1)/num_intermediate_means);            
            if ismember(i_t, floor(im_rng_hi:im_rng_hi:size(var_mtx,1)))
                im_ctr = im_ctr+1;
                
                %timewindow bins in intermediate ranges
                im_range = im_it:i_t;

                %load imean holders
                icell{im_ctr, i_v} = var_mtx(im_range, i_v, :);
                imean(im_ctr, :, i_v) = [mean(var_mtx(im_range, i_v, axes_idx(1))) mean(var_mtx(im_range, i_v, axes_idx(2))) mean(var_mtx(im_range, i_v, axes_idx(3)))];
                icolor_hold(im_ctr, :, i_v) = mean([color_code; im_cc(i_v,:)]);

                %plot imean dot
                plot3(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))), mean(var_mtx(im_range, i_v, axes_idx(3))),...
                    '.', 'Color', icolor_hold(im_ctr, :, i_v), 'MarkerSize', medium_dot_size)
                
                    %encircle imean dot
                    if ismember(variables(i_v), [3 4]) %red for errors
                        plot3(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))), mean(var_mtx(im_range, i_v, axes_idx(3))),...
                            'ro', 'MarkerSize', medium_dot_size/3.5, 'LineWidth', .75)
                    else %black for corrects
                        plot3(mean(var_mtx(im_range, i_v, axes_idx(1))), mean(var_mtx(im_range, i_v, axes_idx(2))), mean(var_mtx(im_range, i_v, axes_idx(3))),...
                            'ko', 'MarkerSize', medium_dot_size/3.5, 'LineWidth', .75)
                    end
                %reset counter
                im_it = i_t+1;
                im_cc(i_v,:) = color_code;
            end
            
            zlabel(['pc' num2str(pc_axes(3))])
        end
    end
    
    %plot means and connect imeans
    if  num_axes ==2
        %plot MEANS 2d
            %{
            plot(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))),...
                '.', 'Color', mean([color_code; base(i_v,:)]), 'MarkerSize', large_dot_size)

                %encircle mean dot
                if ismember(variables(i_v), [3 4]) %red for errors
                    plot(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))),...
                        'ro', 'MarkerSize', large_dot_size/3.5, 'LineWidth', 1.25)
                else %black for corrects
                    plot(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))),...
                        'ko', 'MarkerSize', large_dot_size/3.5, 'LineWidth', 1.25)
                end
            %}
        
        %connect imeans
        for i_im = 1:(size(imean,1)-1)
            plot(imean([i_im i_im+1], 1, i_v), imean([i_im i_im+1], 2, i_v), '-', 'Linewidth', 3, 'Color', mean(icolor_hold([i_im i_im+1], :, i_v)))

            

        end
        
    elseif num_axes==3

        %plot MEANS 3d
        %{
        plot3(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))), mean(var_mtx(:, i_v, axes_idx(3))),...
            '.', 'Color', mean([color_code; base(i_v,:)]), 'MarkerSize', large_dot_size)

            %encircle mean dot
            if ismember(variables(i_v), [3 4]) %red for errors
                plot3(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))), mean(var_mtx(:, i_v, axes_idx(3))),...
                    'ro', 'MarkerSize', large_dot_size/3.5, 'LineWidth', 1.25)
            else %black for corrects
                plot3(mean(var_mtx(:, i_v, axes_idx(1))), mean(var_mtx(:, i_v, axes_idx(2))), mean(var_mtx(:, i_v, axes_idx(3))),...
                    'ko', 'MarkerSize', large_dot_size/3.5, 'LineWidth', 1.25)
            end

            %lable third axis
            zlabel(['pc' num2str(pc_axes(3))])
        %}

        %connect imeans
        for i_im = 1:(size(imean,1)-1)
            plot3(imean([i_im i_im+1], 1, i_v), imean([i_im i_im+1], 2, i_v), imean([i_im i_im+1], 3, i_v), '-', 'Linewidth', 3, 'Color', mean(icolor_hold([i_im i_im+1], :, i_v)))
        end
    
    end
    
end
%label first two axes
xlabel(['pc' num2str(pc_axes(1))])
ylabel(['pc' num2str(pc_axes(2))])

%{
shufs = 1000;
shuffle_mtx = nan(shufs,num_intermediate_means); 
for i_im = 1:num_intermediate_means
    %imean_dists(i_im) = dist(imean(i_im, :, 1), imean(i_im, :, 2)');
    %zc_dists(i_im) = mahal_2cluster_dist(squeeze(icell{i_im, 1}), squeeze(icell{i_im, 2}));
    zc_dists(i_im) = zscore_2cluster_dist(squeeze(icell{i_im, 1}), squeeze(icell{i_im, 2}));
    
    %shuffles
    num_tw = size(squeeze(icell{i_im, 1}),1);
    combo_mtx = [squeeze(icell{i_im, 1}); squeeze(icell{i_im, 2})];
    
    for ishuf = 1:shufs
        
        combo_mtx_idx = randperm(size(combo_mtx,1));
        clustera = combo_mtx(combo_mtx_idx(1:num_tw), :);
        clusterb = combo_mtx(combo_mtx_idx(num_tw+1:end), :);
        
        %shuffle_mtx(ishuf,i_im) = mahal_2cluster_dist(clustera, clusterb);
        shuffle_mtx(ishuf,i_im) = zscore_2cluster_dist(clustera, clusterb);
    end
end

shuffle_mtx = sort(shuffle_mtx);


%figure; bar(imean_dists); title euclidmean
figure; bar(zc_dists); title zdist; 
hold on
plot(1:num_intermediate_means, shuffle_mtx(ceil(.975*shufs),:), 'k-')
%ylim([1 4.1])
%}



