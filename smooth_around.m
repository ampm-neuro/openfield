function vect = smooth_around(vect, smooth_window)
%smooth around the xaxis
%differs from circ smooth in that the y axis is not assumed to vecte circular
%(0:360), only the xaxis is.

    l_vect = length(vect);
    vect = [vect vect vect]; 
    vect = smooth(vect, smooth_window); 
    vect = vect(l_vect+1 : l_vect*2); 

end