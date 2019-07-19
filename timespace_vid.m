function frame = timespace_vid(histc_plot)
%plots decoding (proportion of classifications) over four contexts (BBWW) 
%for all spatial pixels
%
%histc_plot comes from pop_timespace_code

%pixles per context
ppc = size(histc_plot,1)/4;

%pixles on edge of context
pec= sqrt(ppc);

%caxis range
cr = [0 .20];


%plot
h = figure; axis off

reading_order = reshape(1:144,6,144/6)';
ro1 = fliplr(reading_order(1:6,:)); ro1 = ro1(:);
ro2 = fliplr(reading_order(7:12,:)); ro2 = ro2(:);
ro3 = fliplr(reading_order(13:18,:)); ro3 = ro3(:);
ro4 = fliplr(reading_order(19:24,:)); ro4 = ro4(:);
reading_order = [ro1;ro2;ro3;ro4];

%to load getframe
bin_count = 0;

%prep video
v = VideoWriter('timespace.avi');
v.FrameRate = 8;
v.Quality = 100;
open(v)

for bin = reading_order'
bin_count = bin_count+1;


    %b1
    subplot(1,4,1); hold on
    imagesc(reshape(histc_plot(1:ppc, bin), pec, pec))
    axis square
    axis off
    axis([.5 pec+.5 .5 pec+.5])
    caxis(cr)

    %b2
    subplot(1,4,2); hold on
    imagesc(reshape(histc_plot(ppc+1:2*ppc, bin), pec, pec))
    axis square
    axis off
    axis([.5 pec+.5 .5 pec+.5])
    caxis(cr)

    %w1
    subplot(1,4,3); hold on
    imagesc(reshape(histc_plot(2*ppc+1:3*ppc, bin), pec, pec))
    axis square
    axis off
    axis([.5 pec+.5 .5 pec+.5])
    caxis(cr)

    %w2
    subplot(1,4,4); hold on
    imagesc(reshape(histc_plot(3*ppc+1:4*ppc, bin), pec, pec))
    axis square
    axis off
    axis([.5 pec+.5 .5 pec+.5])
    caxis(cr)

    
    
    %add rat position
    switch 1
        case bin >=1 & bin <=ppc
            sp = 1;
        case bin >=ppc+1 & bin <=2*ppc
            sp = 2;
        case bin >=2*ppc+1 & bin <=3*ppc
            sp = 3;
        case bin >=3*ppc+1 & bin <=4*ppc
            sp = 4;
    end
    subplot(1,4,sp)
    
    if rem(bin, ppc) > 0
        [lpi, lpj] = ind2sub([pec pec], rem(bin, ppc));
    else
        lpi = pec;
        lpj = pec;
    end
    plot(lpj, lpi, 'ko', 'markersize', 10)
    
    set(gcf, 'Position', [000, 600, 1000, 400])
    
    %video
    frame(bin_count) = getframe(h);
    
    %gif
    %{
    im = frame2im(frame(bin_count));
    [imind,cm] = rgb2ind(im,256);
    outfile = '/Users/ampm/Documents/MATLAB/PCia/timespace_decode.gif';
    if bin_ct== 1
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
    elseif bin_ct>1
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end
    %}
    
    %video
    writeVideo(v, frame(bin_count));

end

close(v);
end