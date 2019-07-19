function dm = diag_mask(cols)

dm = zeros(cols);

for c = 1:cols
    dm(c,c) = 1;
end
    dm = logical(dm);


end