
file = '1916/R1916_rec12_bwwb_3_21_17'

[eptrials, clusters, context_ids, session_length] = loadnl(file);
%[pos, TimestampsVT, event, flags] = loadnl(file);


%for CSC data, see: 'loadnl'