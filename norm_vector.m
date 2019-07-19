function [normal_vector] = norm_vector(line_p1, line_p2)
%function [plane_out] = perpendicular_plane(line_p1, line_p2)
%
% finds the normal vector connecting points line_p1 and line_p2.
%


%ANSWER projecting data on a hyperplane orthaganol to the vector connecting
%the first and third ITI mean
%ng1d_3 = all_counts_bw_3*null(mean(all_counts_bw_3(ismember(context_bw_3, [12]),:)) - mean(all_counts_bw_3(ismember(context_bw_3, [10]),:)));


%SOLVE FOR NORMAL VECTOR THROUGH POINTS
%
%direction of input line
normal_vector=line_p2-line_p1;

%
%PLOT LINE
% Vector Equation of Line x=P1(1) +dv_line(1)*t;
t=linspace(-1,2);
x=line_p1(1) +normal_vector(1)*t; % X coordinate
y=line_p1(2) +normal_vector(2)*t; % Y coordinate
z=line_p1(3) +normal_vector(3)*t; % Z coordinate
plot3(x,y,z,'LineWidth',2,'Color','r')


%PLOT THE OTHAGANOL PLANE
%
mp = (line_p1 + line_p2) ./ 2;
hold on
t=(0:1:360)';
circle0=[cosd(t) sind(t) zeros(length(t),1)];
r=vrrotvec2mat(vrrotvec([0 0 .0001],normal_vector));
circle=(circle0*r'+repmat(mp,length(circle0),1))./10;
plot3(circle(:,1),circle(:,2),circle(:,3), 'Color', [0.4660    0.6740    0.1880]);
%patch(circle(:,1),circle(:,2),circle(:,3), [0.4660    0.6740    0.1880]);
%}


%uv = data*null(normal_vector);





end
