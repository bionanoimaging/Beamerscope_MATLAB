function [ cpp ] = contrastperpixel( grayImage )
%contrastperpixel Contrast per Pixel Calculation from an input image
%   Detailed explanation goes here
kernel = [-1, -1, -1, -1, 8, -1, -1, -1]/8;
diffImage = conv2(double(grayImage), kernel, 'same');
cpp = mean2(diffImage);

end

