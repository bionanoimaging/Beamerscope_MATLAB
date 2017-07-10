%Program for Image / Picture Quality Measures Calculation

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

%Program Description
%This program is the main entry of the application.
%This program calculates the difference Image/Picture Quality Measures

%Clear Memory & Command Window
clc;
clear all;
close all;

%Read Original & Distorted Images
origImg = imread('H:\Matlab\ImageQualityMeasures\OriginalImages\Lena.bmp');
distImg = imread('H:\Matlab\ImageQualityMeasures\DistortedImages\gaussian_noise_20.bmp');

%If the input image is rgb, convert it to gray image
noOfDim = ndims(origImg);
if(noOfDim == 3)
    origImg = rgb2gray(origImg);
end

noOfDim = ndims(distImg);
if(noOfDim == 3)
    distImg = rgb2gray(distImg);
end

%Size Validation
origSiz = size(origImg);
distSiz = size(distImg);
sizErr = isequal(origSiz, distSiz);
if(sizErr == 0)
    disp('Error: Original Image & Distorted Image should be of same dimensions');
    return;
end

%Mean Square Error 
MSE = MeanSquareError(origImg, distImg);
disp('Mean Square Error = ');
disp(MSE);

%Peak Signal to Noise Ratio 
PSNR = PeakSignaltoNoiseRatio(origImg, distImg);
disp('Peak Signal to Noise Ratio = ');
disp(PSNR);

%Normalized Cross-Correlation 
NK = NormalizedCrossCorrelation(origImg, distImg);
disp('MNormalized Cross-Correlation  = ');
disp(NK);

%Average Difference 
AD = AverageDifference(origImg, distImg);
disp('Average Difference  = ');
disp(AD);

%Structural Content 
SC = StructuralContent(origImg, distImg);
disp('Structural Content  = ');
disp(SC);

%Maximum Difference 
MD = MaximumDifference(origImg, distImg);
disp('Maximum Difference = ');
disp(MD);

%Normalized Absolute Error
NAE = NormalizedAbsoluteError(origImg, distImg);
disp('Normalized Absolute Error = ');
disp(NAE);
