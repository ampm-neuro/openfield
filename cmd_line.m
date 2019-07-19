
time_window = 10; 

%all_counts vars from 
[~, ~, context_bw, context_ia,...
all_counts_bw_z, all_countsz_ia_z] = ...
ALL_countwindow(time_window);

%distances    
method = [1, 2];
tws = (12*60) / time_window;

mtx = all_counts_bw_z;
mtx = mtx(:, ~isnan(mtx(1,:)));
idx = context_bw;
mtx = collapse_time(mtx, idx);
[w_bar, b_bar] = spkct_distances_bw(mtx, idx, [1 2; 3 4], [1 3; 2 4; 2 3; 1 4], method(1), method(2));


mtx = all_countsz_ia_z;
mtx = mtx(:, ~isnan(mtx(1,:)));
idx = context_ia;
mtx = collapse_time(mtx, idx);
[~, ~, a_bar, o_bar, ao_bar] = spkct_distances_ia(mtx, idx, [5; 6; 7; 8], [5 6; 7 8; 5 7; 6 8; 5 8; 6 7], method(1), method(2));

%reshape and average
b_bar = reshape(b_bar, tws, 4); b_bar = mean(b_bar,2); cell_e{1} = b_bar;
w_bar = reshape(w_bar, tws, 4);w_bar = mean(w_bar,2); cell_e{2} = w_bar;
a_bar = reshape(a_bar, tws, 2); a_bar = mean(a_bar,2); cell_e{3} = a_bar;
o_bar = reshape(o_bar, tws, 2); o_bar = mean(o_bar,2); cell_e{4} = o_bar;
ao_bar = reshape(ao_bar, tws, 8); ao_bar = mean(ao_bar,2); cell_e{5} = ao_bar;

%figure
errorbar_barplot(cell_e)
ylim([1 1+(mean(b_bar)-1)*1.15])
xlim([.25 5.75])
xticklabels({'Between', 'Within', 'Arrange', 'Object', 'Arg&Obj'})
xlabel('Context Comparison')
ylabel('Distance')