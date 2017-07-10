function [ tm ] = generateCrossMask( nPixel, R, Mask_L )
%generateCrossMask
%% Generate Mask

[xx yy] = meshgrid(-floor(nPixel/2)-1:floor(nPixel/2)-1);
tm=zeros(length(xx)); % tm= Mask
tm=tm+double(abs(R*yy)<=5*Mask_L)+double(abs(R*xx)<= 5*Mask_L); % Mask_L = 1 ?m
tm=double(tm>1);
tm=tm+double(abs(R*yy)<= Mask_L/2)+double(abs(R*xx)<= Mask_L/2);
tm=double(tm>1);

end

