function newpts = rotate_pts(rotang, pts, center)

% newpts take a list of i j coordinates "pts" and returns them rotated around "center" 
%by the angle "rotang"

% rotang is the rotation angle counter clockwise in degrees
% pts are the points within the matrix that you would like to rotate
% center is the point you would like to rotate the pts about (probably
% comxy)

% newpts is the list of rotated i j coordinates size(2,length(pts))

%translate to radians
rotang = rotang/57.2957795;

%build rotation matrix (some calc II thing)
romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];

%rotation occurs around origin, so we temporarily center all points around the origin
pts(:,1)=pts(:,1)-center(:,1);
pts(:,2)=pts(:,2)-center(:,2);

%apply rotation
newpts=(romat*pts')';

%undo centering
newpts(:,1)=newpts(:,1)+center(:,1);
newpts(:,2)=newpts(:,2)+center(:,2);

end