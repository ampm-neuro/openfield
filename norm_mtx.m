function mtx_n = norm_mtx(mtx)


%zscore
mtx = mtx - repmat(nanmin(mtx), size(mtx,1), 1);
mtx_n = mtx./nanmax(mtx);

end