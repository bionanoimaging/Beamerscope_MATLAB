function frequencymeasure = MaxFrequencyMeasure(Image)

result = log(abs(fftshift(fft2((Image)))));
sx = linspace( -size(result,2)/2, size(result,2)/2, size(result,2));
sy = linspace( -size(result,1)/2, size(result,1)/2, size(result,1));
[XX YY] = meshgrid(sx, sy);
[THETA,RHO] = cart2pol(XX, YY);
RHO = (normminmax(RHO));
RHO1 = im2bw(RHO, 0.1);
RHO2 = im2bw(RHO, 0.6);
A = ((result.*(RHO2-RHO1)));
[M,I] = max(A(:));
[I_row, I_col] = ind2sub(size(A),I);

I_row = I_row - size(Image,1)/2;
I_col = I_col - size(Image,2)/2;

frequencymeasure = sqrt(I_row^2  + I_col^2)