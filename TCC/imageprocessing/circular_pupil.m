function  mask  = circular_pupil( X, Y, xpos, ypos, radius )
%CIRCULAR_PUPIL Create a circle of a given radius at a given position.
%   X - A scaled meshgrid of x positions
%   Y - A scaled meshgrid of y positions
%   xpos - the x position of the center of the circle in real coordinates
%   ypos - the y position of the center of the circle in real coordinates
%   radius - the radius of the circle in real coordinates

mask = double( (X - xpos).^2 + (Y - ypos).^2 <= radius^2 );

end

