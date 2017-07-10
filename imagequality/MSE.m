function MSE= MSE(origImg, distImg);
% http://de.mathworks.com/matlabcentral/answers/231932-is-this-how-to-calculate-mean-square-error-for-two-images
[M, N] = size(origImg);
error = origImg - (distImg);
MSE = sum(sum(error .* error)) / (M * N);
end