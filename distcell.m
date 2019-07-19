function rates = distcell(eptrials, clusters, bin, varargin)
%calculates firing rate of each cell at HD points specified by linespacing
%360degrees by input 'bin.'
%
%context input is CONTEXT ID

dist_col = 12;

%plot by context
if nargin > 100
    contexts = varargin{1};
    if nargin == 5
        dist_col = varargin{2};
    end
    if size(contexts, 1) > size(contexts,2)
        contexts = contexts';
    end
    
    for ctx = contexts
        edges = linspace(0, 1.2, bin);
        
        edges(1) = edges(1)+.01;
        edges(end) = edges(end)-.01;
        
        rates = nan(length(edges)-1, length(clusters)+1);

        for b = 1:length(edges)-1
           spikes = histc(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & eptrials(:,6)==ctx, 4), clusters);
           time = sum(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & eptrials(:,6)==ctx, 4)==1);
           rates(b,1) = mean(edges(b:b+1));
           rates(b,2:end) = spikes./(time*.01);
        end


        %figures
        for i = 2:length(clusters)+1
            figure; 
            plot(rates(:,1), rates(:,i))
            title([num2str(clusters(i-1)) ' ' num2str(ctx)])
            %axis([0 360 0 4*(sum(eptrials(:,4)==clusters(i-1))/eptrials(end,1))])
        end
    end
    
%within context collapse across objects
elseif nargin > 3
contexts = varargin{1};
if size(contexts, 1) > size(contexts,2)
    contexts = contexts';
end

for ctx = contexts
    edges = linspace(0, 1.2, bin);

    edges(1) = edges(1)+.01;
    edges(end) = edges(end)-.01;

    rates = nan(length(edges)-1, length(clusters)+1, 4);
    for dist_col = 10:13
        for b = 1:length(edges)-1
           spikes = histc(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & eptrials(:,6)==ctx, 4), clusters);
           time = sum(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & eptrials(:,6)==ctx, 4)==1);
           rates(b,1, dist_col-9) = mean(edges(b:b+1));
           rates(b,2:end, dist_col-9) = spikes./(time*.01);
        end
    end

    %figures
    for i = 2:length(clusters)+1
        figure; 
        plot(mean(rates(:,1, :),3), mean(rates(:,i, :),3))
        title([num2str(clusters(i-1)) ' ' num2str(ctx)])
        %axis([0 360 0 4*(sum(eptrials(:,4)==clusters(i-1))/eptrials(end,1))])
    end
end
    
    
    
%collapse across contexts     
else

    contexts = 1:4;
    
    edges = linspace(0, 1.2, bin);
    rates = nan(length(edges)-1, length(clusters)+1);

    for b = 1:length(edges)-1
       spikes = histc(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & ismember(eptrials(:,5), contexts), 4), clusters);
           time = sum(eptrials(eptrials(:,dist_col)>=edges(b) & eptrials(:,dist_col)<edges(b+1) & ismember(eptrials(:,5), contexts), 4)==1);
       rates(b,1) = mean(edges(b:b+1));
       rates(b,2:end) = spikes./(time*.01);
    end

    %figures
    for i = 2:length(clusters)+1
        figure; 
        plot(rates(:,1), rates(:,i))
        title([num2str(clusters(i-1)) ' ' num2str(dist_col)])
        %axis([0 360 0 4*(sum(eptrials(:,4)==clusters(i-1))/eptrials(end,1))])
    end


end