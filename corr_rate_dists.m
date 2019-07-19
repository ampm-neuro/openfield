function [withins, alls, object_corrs_alo, obj_appear] = corr_rate_dists(rate_dists)
%computes assorted correlations between rate distributions computed by
%ALL_HD_relative_to_obj
%
%withins is a num_of_cells X num_of_contexts matrix containing the
%similarity (in pearson's r) between firing rate distributions relative to
%the objects in each context. (the average pairwise correlation between HD
%rate distributions to each object).
%
%alls contains the average pairwise correlation between HD rate
%distributions to each object across all contexts. How similar are object
%responses in general.
%
%object_corrs_alo has 3 columns reporting allocentric HD. The first column
%reports basic allocentric HD (relative to north). The second column
%reports the correlation between the two arrangment conditions. The third
%column reports the correlation between the two appearance conditions.
%
%obj_appear col 1 is mean corr between an object in one context and other 3
%objects in same context. col 2 is mean corr between an object in one 
%context and other 3 objects in other context.

%all within-context between-object correlations
withins = nan(length(rate_dists), 4);
%all between-object correlations
alls = nan(length(rate_dists), 1);
%pairwise allocentric hd correlations
object_corrs_alo = nan(length(rate_dists), 3); %all, arrng, appear 
%comparing object responses to other objects in the same or other
%appearance context
obj_appear = nan(length(rate_dists),2); %same, diff ctxt

%for all cells
for cell = 1:length(rate_dists)

    %pairwise object hd correlations
    object_corrs_o_w = nan(length(pairwise(1:4)), 4);
    object_corrs_o_a = nan(length(pairwise(1:16)), 1);

    %pairwise allocentric hd correlations
    object_corrs_alo_local = nan(length(pairwise(1:4)), 1);
    
    %count
    count = 0;
    count_a = 0;
    
    %do object responses vary with appearance?
    obj_appear_comps = nan(4, 2, 8); %all objs (match will be left nan), contexts same/other (7&8), all objs (4objs * 2ctxs)
    
    %iterate through all pairs of FR distribution pairs
    for ctx1 = 1:4
        count_w = 0;
        count_w_alo = 0;
        for ctx2 = 1:4
            for obj1 = 1:5
                for obj2 = 1:5
                    if ~isequal([ctx1 obj1],[ctx2 obj2])
                        count = count+1;
                        
                        %object corrs
                        if ismember(obj1,1:4) && ismember(obj2,1:4)
                            object_corr = corr(...
                                rate_dists{cell}{obj1, ctx1}', rate_dists{cell}{obj2, ctx2}');
                            object_corrs_o_a(count) = object_corr;
                           if ctx1==ctx2
                               count_w = count_w + 1;
                               object_corrs_o_w(count_w, ctx1) = object_corr;
                           else
                               count_w=0;
                           end
                            
                           
                            %do object responses vary with appearance?
                            if ismember(ctx1, 3:4) && ismember(ctx2, 3:4) && ~isequal(obj1,obj2)
                                %all objs (match will be left nan), contexts same/other (7&8), all objs (4objs * 2ctxs)
                                %
                                %each obj has it's own page
                                %rows are pairwise object comparisons
                                %columns are for the two contexts
                                if ctx1==3 
                                    obj_appear_comps(obj2, ctx2-2, obj1) = object_corr;
                                elseif ctx1==4
                                    obj_appear_comps(obj2, setdiff([1 2],ctx2-2), 4+obj1) = object_corr;
                                end
                            end
                           
                        %allocentric corrs
                        elseif obj1==5 && obj2==5
                            count_w_alo = count_w_alo + 1;
                            object_corr = corr(...
                                    rate_dists{cell}{obj1, ctx1}', rate_dists{cell}{obj2, ctx2}');
                                object_corrs_alo_local(count_w_alo) = object_corr;
                            
                            if ctx1==1 && ctx2==2
                                object_corrs_alo(cell, 2) = object_corr;
                            
                            elseif ctx1==3 && ctx2==4
                                object_corrs_alo(cell, 3) = object_corr;
                            
                            end
                        end
                    end
                end
            end
        end
    end
    
    obj_appear_comps
    
obj_appear(cell,:) = nanmean(nanmean(obj_appear_comps), 3);
withins(cell,:) = nanmean(object_corrs_o_w);
alls(cell) = nanmean(object_corrs_o_a);
object_corrs_alo(cell, 1) = nanmean(object_corrs_alo_local);
end
 
end