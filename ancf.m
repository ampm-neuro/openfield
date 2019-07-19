function imp_rows = ancf(eptrials, SPECIFIC_CELL)
%finds rows where eptrials(:,5) numbers change and plot the long context
%visits in a subplot

%preallocate
imp_rows = nan(length(unique(eptrials(:,5)))+1, 1);
imp_rows(end,1) = length(eptrials(:,1));

%counter
counter = unique(eptrials(:,5)); %unique flag numbers
count = 1;

%IMPORTANT ROWS
%for loop to check for number changes
for i1 = 1:length(eptrials(:,5))
    
    %counter(count)
    if count > length(counter) %this is for debugging
        continue
     
    %if we find the next unique flag number
    elseif eptrials(i1,5)==counter(count)
        
        %write down the row number in our preallocated variable
        imp_rows(count,1)=i1;
        %increase counter after finding number change
        count = count+1;
    end
    
end

%PLOTTING WITHIN IMPORTANT ROW RANGES
figure
count=1;
for i2 = 2:length(unique(eptrials(:,5)))+1
    
    %find how long each context visit was
    range_rows = [imp_rows(i2-1,1) imp_rows(i2,1)];
    range = diff(range_rows);
    
    %was context visit long enough to be a "real" context visit (not an
    %ITI)?
    if range > 100000

        %context that passes range test
        context = i2-1;
        
        %plot
        subplot(1,7,count) %subplot determined by counter
        plot(eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 2), eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 3)) %plot only within important row ranges
       
        %ADD SPIKE FIRING OVERLAID
        hold on
        
        %downsample
        how_many_spikes_percent = 5;
        
        idx = eptrials(:,4)==SPECIFIC_CELL & eptrials(:,5)==context;
        idx = idx.*rand(size(idx))>((100-how_many_spikes_percent)/100);
        
        plot(eptrials(idx,2),eptrials(idx,3),'r.')
        hold off
        
        axis([min(eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 2)-20) max(eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 2)+20) min(eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 3)-20) max(eptrials(imp_rows(i2-1,1):imp_rows(i2,1), 3)+20)]) %make great axes
        title(['Context ', num2str(count)],'fontsize', 16) %title here works magic with strings
        
        count = count+1; %increase counter after successful plot
    end
end
end