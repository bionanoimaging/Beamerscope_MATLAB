function [ Ic ] = generateSegments( circle, coefficient, NAo, NAc )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



PixelDiameterNA = sum(circle(round(size(circle,1)/2), :))/2;
pixelsize = NAc/PixelDiameterNA;

dimy = linspace(-pixelsize.*size(circle,1)/2, pixelsize.*size(circle,1)/2, size(circle,1));
dimx = linspace(-pixelsize.*size(circle,2)/2, pixelsize.*size(circle,2)/2, size(circle,2));

[mm, nn]=meshgrid(dimx, dimy);

% Obtain azimuthal and radial co-ordinates in the pupil plane.
%---------------------
[mtheta, mr]=cart2pol(mm,nn);



S=NAc/NAo;
n_segment = 12;      % number of colour segments
n_subsegments = round(size(coefficient,2)/n_segment,0);
Ic = zeros(size(circle));
iter = 1;

for j = 0:n_subsegments-1
    
    for i=0:n_segment-1
        
        Isegment = zeros(size(circle));
        
        % i-th segment in one of the annuli
        i_segment = rem(i, n_segment);
        
        
        
        
        NAc_i = j*(1/n_subsegments);
        NAc_o = (j+1)*(1/n_subsegments);
        
        NAc_i = NAc_i * NAc;
        NAc_o = NAc_o * NAc;
        
        % Creating the annullar shape 0,1,2
        Isegment=(mr>=NAc_i/NAo ... Inner radius.
            & mr<=NAc_o/NAo ... Outer radius.
            );
        
        % scale rotational symmetric ramp 0..1
        segment_area = ((mtheta)/max(max(mtheta)) * round(n_segment/2)) + round(n_segment/2);
        
        % each segment results from the threshold of the grayvalues
        % filtered by the annular shape of the illumination sector
        % 0,1,2
        
        % this is due to aliasing of the pixelated source, otherwise
        % there will be a gap in the source shape
        if(i_segment == n_segment-1)
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1)*1.00001);
        else
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1));
        end
        
        
        % get i-th segment and sum it up; weigh it with coefficient
        segment_area = segment_area.*Isegment;
        
        Isegment = segment_area.*coefficient(iter);
        
        Ic = Ic + Isegment;
        
        iter = iter + 1;
    end
    
end



end

