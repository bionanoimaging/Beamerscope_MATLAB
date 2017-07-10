%Program for Normalized Absolute Error Calculation

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

function NAE = NormalizedAbsoluteError(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

error = origImg - distImg;

NAE = sum(sum(abs(error))) / sum(sum(origImg));