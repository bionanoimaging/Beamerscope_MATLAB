%% Example I demontsrate various shift fucntion with integer shift
close all;
img=imread('new-york.jpg');
figure; 
imshow(img); title('Original image');

shitVec=[20.78, 18.09];
hFig=figure;
set(hFig, 'Name', 'Integer shift- sharp images');
subplot_tight(2,2,1); imshow( circshift(img, round(shitVec)) ); 
title( sprintf('"circshift" Shifted image- repeating bounds.\n Note: statue legs go up the sky.') );
subplot_tight(2,2,2); imshow( imshift(round(shitVec), img) );
title( sprintf('"imshift" Shifted image- zero bounds.\n Note: black margins on top-left') );
subplot_tight(2,2,3); imshow( normshift(img, round(shitVec)) ); 
title( sprintf('"normshift" Shifted image- zero bounds.\n Note: black margins on top-left') );
subplot_tight(2,2,4); imshow( imshift(round(shitVec), img, 'symmetric') ); 
title( sprintf('"imshift" Shifted image- mirror-reflecting bounds.\n Note: sky added, and torch flame multiplied') );

%% Example II demontsrate various shift fucntion with fractional shift
hFig=figure;
set(hFig, 'Name', 'Fractional shift- blurry images');
subplot_tight(2,2,1); imshow( floatingCircShift(img, shitVec) ); 
title( sprintf('"circshift" Shifted image- repeating bounds.\n Note: statue legs go up the sky.') );
subplot_tight(2,2,2); imshow( imshift(shitVec, img) );
title( sprintf('"imshift" Shifted image- zero bounds.\n Note: black margins on top-left') );
subplot_tight(2,2,3); imshow( normshift(img, shitVec) ); 
title( sprintf('"normshift" Shifted image- zero bounds.\n Note: black margins on top-left') );
subplot_tight(2,2,4); imshow( imshift(shitVec, img, 'symmetric') ); 
title( sprintf('"imshift" Shifted image- mirror-reflecting bounds.\n Note: sky added, and torch flame multiplied') );

% Example III- using fractional shift to generate image interpolation
img=imread('peppers.png');
img=imcrop(img, [190, 80, 60, 60]);
img1_1=img;
img1_2=imshift([0, -0.5], img); % Note the negative sign, due to circhift convention
img2_1=imshift([-0.5, 0], img);
img2_2=imshift([-0.5, -0.5], img);

imdDims=size(img);
% nearest neighbour- by hand
imgNewDims=[2*imdDims(1), 2*imdDims(2), imdDims(3)];
imgNN=zeros( imgNewDims, class(img) );
imgNN(1:2:end, 1:2:end, :)=img;
imgNN(1:2:end, 2:2:end, :)=img;
imgNN(2:2:end, 1:2:end, :)=img;
imgNN(2:2:end, 2:2:end, :)=img;

% Linear interpolation neighbour- by hand
imgInterp=zeros( imgNewDims, class(img) );
imgInterp(1:2:end, 1:2:end, :)=img1_1;
imgInterp(1:2:end, 2:2:end, :)=img1_2;
imgInterp(2:2:end, 1:2:end, :)=img2_1;
imgInterp(2:2:end, 2:2:end, :)=img2_2;


figure;
subplot_tight(2,2,1)
imshow(imgNN); title('Nearest neighbour- by simple repeating image');
subplot_tight(2,2,2)
imshow(imgInterp); title('Bilinear interpolation using "imshift"');
subplot_tight(2,2,3)
imshow( imresize(img, 2, 'nearest') ); title('Nearest neighbour interpolation using "imresize"');
subplot_tight(2,2,4)
imshow( imresize(img, 2, 'bilinear') ); title('Bilinear interpolation  using "imresize"');