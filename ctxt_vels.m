function [velocs, cumdists, velocs_first, cumdists_first, velocs_second, cumdists_second] = ctxt_vels(eptrials, contexts)
%function ctxt_rates(eptrials, contexts)
%
%calculates the average velocity and total distance traveled in each of
%the context conditions in contexts.
%
%also calculates the firing rate for first and second halves in each
%context

%preallocate
velocs = nan(1, length(contexts));
cumdists = nan(1, length(contexts));
velocs_first = nan(1, length(contexts));
cumdists_first = nan(1, length(contexts));
velocs_second = nan(1, length(contexts));
cumdists_second = nan(1, length(contexts));


%counters
context_num = 0;

%orientation check
if size(contexts,1) > size(contexts,2)
    contexts = contexts';
end
    
for trial_type = contexts

    %count contexts
    context_num = context_num + 1;

    %ALL
        %distance
        dist_hold = eptrials(eptrials(:,6)==trial_type & eptrials(:,4)==1,2:3);
        cumdistance =[0;cumsum(sqrt(diff(dist_hold(:,1)).^2 + diff(dist_hold(:,2)).^2))];
        cumdist = cumdistance(end);
        
        %veloc 7, acceleration 8
        %veloc = cumdist/(length(cumdistance)*.01);
        veloc = nanmean(abs(eptrials(eptrials(:,6)==trial_type,8)));
       
        %load
        velocs(context_num) = veloc;
        cumdists(context_num) = cumdist;

    %FIRST
%{
        %halfway point
        hp = mean([min(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type, 1)) max(eptrials(eptrials(:,4)==1 & eptrials(:,6)==trial_type, 1))]);

        %distance
        dist_hold = eptrials(eptrials(:,5)==trial_type & eptrials(:,4)==1 & eptrials(:,1) < hp,2:3);
        cumdistance =[0;cumsum(sqrt(diff(dist_hold(:,1)).^2 + diff(dist_hold(:,2)).^2))];
        cumdist = cumdistance(end);
        
        %veloc
        veloc = cumdist/(length(cumdistance)*.01);
       
        %load
        velocs_first(context_num) = veloc;
        cumdists_first(context_num) = cumdist;


    %SECOND
        %distance
        dist_hold = eptrials(eptrials(:,6)==trial_type & eptrials(:,4)==1 & eptrials(:,1) > hp,2:3);
        cumdistance =[0;cumsum(sqrt(diff(dist_hold(:,1)).^2 + diff(dist_hold(:,2)).^2))];
        cumdist = cumdistance(end);
        
        %veloc
        veloc = cumdist/(length(cumdistance)*.01);
       
        %load
        velocs_second(context_num) = veloc;
        cumdists_second(context_num) = cumdist;
%}
end

end