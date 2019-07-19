function test_hd(eptrials, speed, varargin)

%test_hd takes vectors of x positions, y positions, and head direction, and
%plays a video of the rat moving around the maze.
%
%eptrials is only used to make the maze boundaries, and can be deleted.
%
%quiver_arrow_size is only asthetic and can be deleted. its based on 3rd
%party code. google: adjust_quiver_arrowhead_size

%if input modified column (9 is HD)
if nargin == 3
    col = varargin{1};
else
    col = 9;%use HD column
end


figure;
object_grid(mode(eptrials(:,6)))
axis([-0.005 1.005 -0.005 1.005]); 
axis square
set(gca,'Xtick',[]);set(gca,'Ytick',[]);
set(gca,'TickLength',[0, 0]);
box on


%sections(eptrials);
axis([-0.005 1.005 -0.005 1.005])

%iterate through sample time points
for pos = 1:speed:size(eptrials,1) %1:length(sample(:,1))
    
    %magic arrow ingredients (see help of quiver)
    angle = eptrials(pos,col);
    s = eptrials(pos,2);
    t = eptrials(pos,3);
    u = sin(angle*(pi/180));
    v = cos(angle*(pi/180));
    line_length = .1;
    arrow_size = 4;    
    
    %plot pos and hd
    h1 = quiver(s, t, u, v, line_length, 'k', 'linewidth', 3, 'Marker', 'o', 'MaxHeadSize', 1000, 'AutoScale', 'off');

    %time legend
    frame = getframe;
    frame2im(frame);
    delete(h1);
    %delete(h2);
    
end


end