function [I] = calculateImage(self, object)
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

disp('%%%%%%Calculate Aerial Image%%%%%%%')
tic
objectspectrum=fftshift(ifft2(ifftshift(object)));
I = zeros(size(objectspectrum));
%padsize = round((size(objectspectrum, 2) - size(self.eigenfunction(:, :, 1),2))/2);
scalefactor = size(I,2)/size(self.eigenfunction(:, :, 1),2);

% Finally aerial image can be calculated by equation (3):
for i=1:self.N
    %if(i<=size(self.eigenvalue,2))
        
        %case 0 - padd pupil..wrong! because fft scales its input, physical
        %relationship between spatial and frequency coordinates is not
        %ensured anymore!
        % eigfun = padarray(self.eigenfunction(:, :, i), [ padsize, padsize]);

        
        %case 1 interpolate (upscale) pupil
        eigfun = imresize(self.eigenfunction(:, :, i), scalefactor);
        
        aerial=objectspectrum.*eigfun;%self.eigenfunction(:, :, i) ;
        FTaerial=fftshift(fft2(ifftshift(aerial)));
        I=I+self.eigenvalue(i).^2*abs(FTaerial).^2;
        %         imagesc(log(abs(self.eigenfunction(:, :, i))))
        %         axis square
   % end    
end



toc
