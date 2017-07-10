function [ NormError ] = NRMSE( I1, I2 )
%NRMSE Calculates the normalized Root Mean Squared Error 
%   I1, I2 are two images which needs to be compared
    NormError =  sqrt(mean(((I1(:) - I2(:)).^2)))*100;

end

