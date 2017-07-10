function  mask  = semicircular_pupil( X, Y, xpos, ypos, radius, angle )
%CIRCULAR_PUPIL Create a circle of a given radius at a given position.
%   X - A scaled meshgrid of x positions
%   Y - A scaled meshgrid of y positions
%   xpos - the x position of the center of the circle in real coordinates
%   ypos - the y position of the center of the circle in real coordinates
%   radius - the radius of the circle in real coordinates
%   angle - rotation angle of sem circle



    mask = ((X - xpos).^2 + (Y - ypos).^2 < radius.^2) & ((Y - ypos)*cos(angle)-(X - xpos)*sin(angle)>=0); 

end

