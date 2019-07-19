function pairs_out = pairwise(vect_in)
%outputs each unique pairwise pairing (order does not matter)

vect_in = unique(vect_in);

perms = (factorial(length(vect_in)))...
    / factorial((length(vect_in) - 2));
perms = round(perms);

pairs_out = nan(perms, 2);

count_row_idx = 0;
for ivectin = 1:length(vect_in)
    
    not_ivectin = setdiff(vect_in, vect_in(ivectin)); 
    
    for iitem = 1:length(not_ivectin)
    
        count_row_idx = count_row_idx+1;
        
        pairs_out(count_row_idx, 1) = vect_in(ivectin);
        pairs_out(count_row_idx, 2) = not_ivectin(iitem);
        
    end
end