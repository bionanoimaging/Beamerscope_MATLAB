function I = calculateImageWithZernikePoly(x, objectspectrum, eigenfunction, eigenvalue, orderKernel)

I = zeros(size(objectspectrum));
% Finally aerial image can be calculated by equation (3):

for(j = 1:size(x,2))
for i=1:orderKernel
aerial=objectspectrum.*eigenfunction(:, :, i, j);
FTaerial=fftshift(fft2(ifftshift(aerial)));
I=I+x(j)*eigenvalue(i, j).^2*abs(FTaerial).^2;
end
end

end



