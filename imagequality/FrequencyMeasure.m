function frequencymeasure = FrequencyMeasure(Image)

result = log(abs(fftshift(fft2((Image)))));
sx = linspace( -size(result,2)/2, size(result,2)/2, size(result,2));
sy = linspace( -size(result,1)/2, size(result,1)/2, size(result,1));
[XX YY] = meshgrid(sx, sy);
[THETA,RHO] = cart2pol(XX, YY);
RHO = (normminmax(RHO));
RHO = im2bw(RHO, 0.2);
frequencymeasure = sum(sum(result.*RHO));
