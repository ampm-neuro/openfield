function videoread_HD(eptrials, video_file, sample_rate)
%makes a video of the four context visits overlayed with a head direction
%arrow.


    % Access video file
    v = VideoReader(video_file);

    % Preallocate structure to store video frames
    % Note that this is just an initial guess for preallocation based on
    % duration and framerate, the video may have fewer or more frames
    nFrames = ceil(v.FrameRate*v.Duration);
    s(nFrames) = struct('cdata',[],'colormap',[]);

    % Set up figure, axes to hold image, image, and plot
    % This allows us to keep everything positioned in the proper order, and
    % just change the image data and the plot location every loop iteration
    hFig = figure('MenuBar','none',...
        'Units','pixels',...
        'Position',[100 100 v.Width v.Height]);
    hAx = axes('Parent',hFig,...
        'Units','pixels',...
        'Position',[0 0 v.Width v.Height],...
        'NextPlot','add',...
        'Visible','off',...
        'XTick',[],...
        'YTick',[]);
    hIm = image(uint8(zeros(v.Height,v.Width,3)),...
        'Parent',hAx);
    hLine(1) = plot(hAx,[1 v.Width],[3 3],'-b','LineWidth',2);
    hLine(2) = plot(hAx,[1 1],[3 3],'-b','LineWidth',4);
    hLine(3) = plot(hAx,1,3,'ob','MarkerSize',10,'MarkerFaceColor','b');

    % Loop through video, grabbing frames and updating plots
    k = 1;
    while hasFrame(v)
        im = readFrame(v);
        hIm.CData = im;
        % Simulate progress bar
        hLine(2).XData(2) = v.Width*k/nFrames;
        hLine(3).XData = v.Width*k/nFrames;
        drawnow
        % Save the frame in structure for later saving to video file
        s(k) = getframe(hAx);
        k = k+1;
    end

    % Remove any unused structure array elements
    s(k:end) = [];

    % Open a new figure and play the movie from the structure
    hFig2 = figure;
    movie(hFig2,s,1,v.FrameRate);

    % Write to a video file
    % This could be done within the original loop, but I wanted to show it
    % separately
    vOut = VideoWriter('xylophoneProgressBar.mp4','MPEG-4');
    vOut.FrameRate = v.FrameRate;
    open(vOut)
    for k = 1:numel(s)
        writeVideo(vOut,s(k))
    end
    close(vOut)


end