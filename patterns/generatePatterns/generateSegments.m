function [ Ic ] = generateSegments( sizeMask, n_segments, n_rings )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


dimy = linspace(-sizeMask(1)/2, sizeMask(1)/2, sizeMask(2));
dimx = linspace(-sizeMask(1)/2, sizeMask(1)/2, sizeMask(2));

[mm, nn]=meshgrid(dimx, dimy);

% Obtain azimuthal and radial co-ordinates in the pupil plane.
%---------------------
[mtheta, mr]=cart2pol(mm,nn);

mr = mr./max(max(mr));


Ic = zeros(sizeMask(1), sizeMask(2), n_segments*n_rings);
iter = 1;

for j = 0:n_rings-1
    
    for i=0:n_segments-1
        
        % create segment layer
        Isegment = zeros(sizeMask);
        
        % i-th segment in one of the annuli
        i_segment = rem(i, n_segments);
        
        NAc_i = (j*(1/n_rings)).^1.5;
        NAc_o = ((j+1)*(1/n_rings)).^1.5;
        
        % Creating the annullar shape 0,1,2
        Isegment=(mr>=NAc_i ... Inner radius.
            & mr<=NAc_o ... Outer radius.
            );

        % scale rotational symmetric ramp 0..1
        segment_area = ((mtheta)/max(max(mtheta)) * round(n_segments/2)) + round(n_segments/2);
        
        % each segment results from the threshold of the grayvalues
        % filtered by the annular shape of the illumination sector
        % 0,1,2
        
        % this is due to aliasing of the pixelated source, otherwise
        % there will be a gap in the source shape
        if(i_segment == n_segments-1)
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1)*1.00001);
        else
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1));
        end
        
        % get i-th segment and sum it up; weigh it with coefficient
        Ic(:,:,iter)  = segment_area.*Isegment;
        
        iter = iter + 1;
    end
end

end

