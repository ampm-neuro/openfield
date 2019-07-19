function object_loc = object_grid(context, varargin)
%plot objec tlocations in PCIA


if nargin == 2
    bins = varargin{1};
    figure_on = 1;
elseif nargin == 3
    bins = varargin{1};
    figure_on = varargin{2};
else
    bins = 1;
    figure_on = 1;
end

%object radius
radius = .17/2;
radius = radius*bins;

if figure_on == 1
    hold on
    axis square
    if bins == 1;
        axis([0 1 0 1])
    else
        axis([.5 bins+.5 .5 bins+.5])
    end
end
%buffer
%{
plot([0 1], [0.125 0.125], 'k-')
plot([0 1], [1-0.125 1-0.125], 'k-')
plot([0.125 0.125], [0 1], 'k-')
plot([1-0.125 1-0.125], [0 1], 'k-')

%grid
plot([0 1], [0.3125 0.3125], 'k-')
plot([0 1], [0.5 0.5], 'k-')
plot([0 1], [0.6875 0.6875], 'k-')
plot([0.3125 0.3125], [0 1], 'k-')
plot([0.5 0.5], [0 1], 'k-')
plot([0.6875 0.6875], [0 1], 'k-')
%}


%define object centers
if context == 5
    object_loc(1, :) = [mean([0.5 0.6875]) mean([0.6875 0.875])] + [0.015 0.02];
    object_loc(2, :) = [mean([0.3125 0.5]) mean([0.5 0.6875])] + [0.03 0.02];
    object_loc(3, :) = [mean([0.6875 1-0.125]) mean([0.3125 0.5])] + [0.01 0.04];
    object_loc(4, :) = [mean([0.3125 0.5]) mean([0.125 0.3125])] + [0.03 0.02];
    
elseif context == 6
    object_loc(1, :) = [mean([0.5 0.6875]) mean([0.5 0.6875])] + [-0.0025 0.03];
    object_loc(2, :) = [mean([0.6875 1-0.125]) mean([0.5 0.6875])] + [0 0.03];
    object_loc(3, :) = [mean([0.125 0.3125]) mean([0.3125 0.5])] + [0.02 0.03];
    object_loc(4, :) = [mean([0.3125 0.5]) mean([0.3125 0.5])] + [0.02 0.03];
    
elseif ismember(context,[7 8])
    object_loc(1, :) = [mean([0.125 0.3125]) mean([0.6875 1-0.125])] + [0.02 0.01];
    object_loc(2, :) = [mean([0.6875 1-0.125]) mean([0.6875 1-0.125])] + [0.02 0.01];
    object_loc(3, :) = [mean([0.125 0.3125]) mean([0.125 0.3125])] + [0.02 0.02];
    object_loc(4, :) = [mean([0.6875 1-0.125]) mean([0.125 0.3125])] + [0.02 0.02];
    
else
    error('context not specified')
    
end
object_loc = object_loc.*bins;
    
%report centers
%object_loc

figure; colors = get(gca,'ColorOrder'); close
%colors = zeros(size(colors));

%plot centers and bounds of each object
if figure_on==1
    for i = 1:4

        circle(object_loc(i,:), radius, 3, colors, i)

        rect_pos = [object_loc(i,1)-radius object_loc(i,2)-radius radius*2 radius*2];
        rectangle('Position', rect_pos, 'Curvature', [1 1], 'FaceColor', [1 1 1])

        %custom objects
        if ismember(context, [5 6])
            %plot(object_loc(i,1), object_loc(i,2), 'k.', 'Markersize', 50)
            plot(object_loc(i,1), object_loc(i,2), 'color', colors(i,:), 'Markersize', 40)
            circle(object_loc(i,:), radius*0.65, 4, colors, i)
        elseif context == 7
            %plot(object_loc(i,1), object_loc(i,2), 'k.', 'Markersize', 50)
            plot(object_loc(i,1), object_loc(i,2), 'color', colors(i,:), 'Markersize', 40)
            circle(object_loc(i,:), radius*0.75, 4, colors, i)
            %circle(object_loc(i,:), radius*0.475, 4)
        elseif context == 8
            %plot(object_loc(i,1), object_loc(i,2), 'k.', 'Markersize', 100)
            plot(object_loc(i,1), object_loc(i,2), 'color', colors(i,:), 'Markersize', 40)
            plot([object_loc(i,1) object_loc(i,1)], [object_loc(i,2)-radius object_loc(i,2)+radius], 'color', colors(i,:), 'linewidth', 4)
            plot([object_loc(i,1)-radius object_loc(i,1)+radius], [object_loc(i,2) object_loc(i,2)], 'color', colors(i,:), 'linewidth', 4)
        end

    end
end

function circle(xy,r,linewidth, colors, i)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
    hold on
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    %plot(xy(1)+xp,xy(2)+yp, 'k', 'linewidth', linewidth);
    plot(xy(1)+xp,xy(2)+yp, 'color', colors(i,:), 'linewidth', linewidth);

end 
    
end
