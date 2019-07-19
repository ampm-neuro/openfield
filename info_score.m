function scores = info_score(eptrials, bins, varargin)
%finds the spatial information content of all clusters in bits per spike

if nargin ==3
    clusters = varargin{1};
else
    clusters = unique(eptrials(~isnan(eptrials(:,4)),4))';
end

%how much time the rat spent in every pixle
min_time = 0.2; %seconds

%preallocate
scores = nan(size(clusters));
info_content = nan(bins);

for clust = clusters
    
    %mean spike counts and occupancy (in seconds) for every bin
    [~, spk_ct, spc_occ] = ...
        trlfree_heatmap(eptrials, clust, bins, 0);
    
    %set non-visted pixles to nans for both spikes and occupancy
    spc_occ(spc_occ<min_time) = 0;
    spk_ct(spc_occ<min_time) = 0;
    spc_occ(spc_occ==0) = nan;
    spk_ct(spk_ct==0) = nan;
    
    %check for number of visited pixles
    %prop_pix_vis = sum(sum(~isnan(spc_occ)))/bins^2
    %if prop_pix_vis < .09*bins^2   ||  prop_pix_vis > .22*bins^2 
    %    continue
    %end
        
    %spatial probabilities
    spc_p = spc_occ./nansum(spc_occ(:));
    
    %rate heatmap
    %spk_rate = spk_ct./spc_occ;
    %spk_rate = smoothmtx(spk_rate, 7);
    %
    spk_rate = skagg_smooth(spk_ct, spc_occ);
    %spk_rate = smooth2a(spk_rate, 1);

    %R is the overall mean firing rate
    R = nansum(nansum(spk_rate.*spc_p));
    
    %exclude cells firing less than 3hz
    if R < 3
        continue
    end
    
    
    
    
    %heatmap figure
    %{
    fig_counter = 0;
    if stage==4 && rand(1) > .98
        figure; imagesc(spk_rate);colormap jet; colorbar; caxis([0 R*2]); title rate
        %figure; imagesc(spc_p);colormap jet; colorbar; title spaceprob
        fig_counter = 1;
    elseif stage<4 && rand(1) > .98
        figure; imagesc(spk_rate);colormap jet; colorbar; caxis([0 R*2]); title rate
        %figure; imagesc(spc_p);colormap jet; colorbar; title spaceprob
        fig_counter = 1;
    end
    %}
    
    %calculate info content from rate map
    for i = 1:numel(spc_occ)
        Pi = spc_p(i);%probability of occupancy of bin i
        Ri = spk_rate(i);%mean firing rate for bin i
        rr = Ri/R;%relative rate
        info_content(i) = Pi*(rr)*log2(rr);%load info content for each bin
    end
    
    
    
    %info content figure
    %{
    if fig_counter ==1
        figure; imagesc(info_content);colormap jet; colorbar; caxis([-.05 .05])
        title(['infocontent ' num2str(nansum(info_content(:))) ' ' num2str(stage)])
    end
    %}
    
    %load cluster info scores
    scores(clusters==clust) = nansum(nansum(info_content)); 
    
    
     %figure; imagesc(spk_rate);colormap jet; colorbar; caxis([0 R*2]); title(num2str(scores(clusters==clust)))
    
        
    %clear info_content for next cluster
    info_content = nan(bins);
    
    
   
    
end


%smooth function
    function smtx_out = smoothmtx(mtx_in, smooth_factor)
        %smtx_out = mtx_in;
        size_fix = mtx_in; 
        size_fix(~isnan(size_fix))=1;

        mask = fspecial('Gaussian',[smooth_factor smooth_factor],1.5);
        
        smtx_out = conv2nan(mtx_in, mask, 'same');smtx_out = smtx_out.*size_fix;
        %smtx_out = conv2nan(smtx_out, mask, 'same');smtx_out = smtx_out.*size_fix;
        %smtx_out = conv2nan(smtx_out, mask, 'same');smtx_out = smtx_out.*size_fix;
    end

end