function [ output_args ] = time_veloc( mtx, idx )
%plots change in rate / time for consequtive time windows over all sessions
%   plots time code (drift) as changing velocity over the session. 
%
%   Because ITIs and session visits are inequal between rats, this 
%   analysis takes the 30s before and 6min after entrance to each context 
%   6min before and 30s after exit from each context
%

%context order
ctx_ids = unique(idx);
ctx_odr = ctx_ids([5 1 6 2 7 3 8 4 9]);






end

