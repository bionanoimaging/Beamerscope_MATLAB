%Program for Peak Signal to Noise Ratio Calculation

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

function PSNR = PeakSignaltoNoiseRatio(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;
MSE = sum(sum(error .* error)) / (M * N);

if(MSE > 0)
    PSNR = 10*log(255*255/MSE) / log(10);
else
    PSNR = 99;
end