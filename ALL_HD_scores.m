function [mis_all_allo, mis_all_obj, rays_all_allo, rays_all_obj, ...
    rs_all_allo, rs_all_obj, obj_hd_id, allo_hd_id] = ALL_HD_scores(rd, allo_hd)
%calculate miller info scores and rayleigh scores for every neuron
%
%try loading 'HD_bw_distribs.mat'

%ALL OBJECT CONDITIONS
mis_all_obj = nan(length(rd),1); 
mis_all_allo = nan(length(rd),1); 
rays_all_obj = nan(length(rd),1);
rays_all_allo = nan(length(rd),1); 
rs_all_obj = nan(length(rd),2);
rs_all_allo = nan(length(rd),2);
obj_hd_id = nan(length(rd),1);
allo_hd_id = nan(length(rd),1);


%colors
colors_woo = [...
0    0.4470    0.7410;... %blue for arrange 1
0.9290    0.6940    0.1250;... %purple for arrange 2
0.8500    0.3250    0.0980;... %red-orange for appear 1
0.4940    0.1840    0.5560]; %yellow for appear 2


%for all cells in rd
for coi = 1:length(rd)
    ds_obj = [];
    ds_allo = [];
    
    %cols and rows of rd{}
    for i1 = 1:5
        for i2 = 1:4
            if i1~=5
                %all obj cols/rows
                ds_obj = [ds_obj; rd{coi}{i2,i1}];
            end
        end
    end

    %obj
    mis_all_obj(coi) = miller_infoscore(ds_obj);
    rs_all_obj(coi, :) = pairwise_corrs(ds_obj);
    %[~, ray_R] = rayleigh(nanmean(ds_obj));
    %rays_all_obj(coi) = ray_R;
    
    %allo
    ahd = nan(4, size(allo_hd,2));
    for i = 1:4
        ahd(i,:) = allo_hd(coi,:,i);
    end
    %ahd = correct_means(ahd);
    mis_all_allo(coi) = miller_infoscore(ahd);
    %[~, ray_R] = rayleigh(nanmean(ds_allo));
    rs_all_allo(coi, :) = pairwise_corrs(ahd);
    %rays_all_allo(coi) = ray_R;
    
    
    %examples of outliers
    %{
    if mis_all_allo(coi)<0.1 && mis_all_obj(coi)<0.1
        figure;
        for ctx_id = 1:4
            subplot(1,2,1)
            polarplot(ds_obj((4*(ctx_id-1)+1):(4*(ctx_id-1)+1)+3,:)', 'color', colors_woo(ctx_id,:)) 
            title([num2str(coi) ' obj ' num2str(mis_all_obj(coi))])
            hold on;
        end
         
        for ctx_id = 1:4
            subplot(1,2,2)
            polarplot(ahd(ctx_id,:), 'color', colors_woo(ctx_id,:))
            title([num2str(coi) ' allo ' num2str(mis_all_allo(coi))])
            hold on; 
        end
    end
    %}
    
    %shuffles
    [hishuf_obj, loshuf_obj] = ...
    miller_infoscore_shuffle(ds_obj, 10, .05);

    if mis_all_obj(coi)>hishuf_obj
        obj_hd_id(coi) = 1;
    elseif mis_all_obj(coi)<loshuf_obj
        obj_hd_id(coi) = -1;
    else
        obj_hd_id(coi) = 0;
    end


    [hishuf_allo, loshuf_allo] = ...
    miller_infoscore_shuffle(ahd, 10, .05);

%figure; plot(distr)
    if mis_all_allo(coi)>hishuf_allo
        allo_hd_id(coi) = 1;
    elseif mis_all_allo(coi)<loshuf_allo
        allo_hd_id(coi) = -1;
    else
        allo_hd_id(coi) = 0;
    end
end
    figure; plot(mis_all_allo, mis_all_obj, 'o')
end


function ds = correct_means(ds)
%set all means equal to highest mean
    dsmeans = mean(ds,2);
    himean = max(dsmeans);
    meandiffs = abs(dsmeans-himean);
    ds = ds+meandiffs;
end

function out = pairwise_corrs(ds)
out = [];
for ids1 = 1:size(ds,1)
    for ids2 = 1:size(ds,1)
        if ids1~=ids2
            [corrR, corrP] = corr(ds(ids1,:)', ds(ids2,:)');
            out = [out; [corrR, corrP]];
        end
    end
end
out = [mean(out(:,1)) max(out(:,2))];

end