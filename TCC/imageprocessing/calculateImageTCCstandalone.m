function [I] = calculateImageTCCstandalone(eigenfunction, eigenvalue, object, M, N)
% This function calculates the aerial image from the objects spectrum and 
% the TCC according to the paper:
% Aerial Image Simulation for partial coherent system with programming 
% development in MATLAB
%
% 
% input:
% - TCC - 2D TCC representation; each TCC is stacked into one row
% - eigenfunction
% - eigenvalue
% - objectspectrum - spectrum of object calculated by using FFT
% - M, N - Dimensions of the TCC
%
%
% output:
% - I - simulated image


%disp('%%%%%%Calculate Aerial Image%%%%%%%')
objectspectrum=fftshift(ifft2(ifftshift(object)));
tic
I = zeros(size(objectspectrum));
% Finally aerial image can be calculated by equation (3):
for i=1:N
aerial=objectspectrum.*eigenfunction(:, :, i);
FTaerial=fftshift(fft2(ifftshift(aerial)));
I=I+eigenvalue(i, i).^2*abs(FTaerial).^2;
end
toc

end
