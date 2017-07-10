%Program for Average Difference Calculation

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

function AD = AverageDifference(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;

AD = sum(sum(error)) / (M * N);