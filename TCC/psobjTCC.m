function y = psobjTCC(x, xmask, objectspectrum, eigenfunction, eigenvalue, orderKernel, object, qualitymethod)
%% Create objective function for the optimization process using the eigenvalues/eigenfunctions created with TCC class
% 
% 
% The software is for private use only and gives no guarrantee, that it's
% working as it should! 
% 
%
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.
%% inputs
% x
% xmask
% objectspectrum
% eigenfunction
% eigenvalue
% orderKernel
% object
% qualitymethod

I = zeros(size(objectspectrum));
% Finally aerial image can be calculated by equation (3):
% padsize = round((size(objectspectrum, 2) - size(eigenfunction(:, :, 1),2))/2);
scalefactor = size(I,2)/size(eigenfunction(:, :, 1),2);


xindex=find(xmask);

for(j = 1:find(x,2))
    for i=1:orderKernel
        %     eigfun = padarray(eigenfunction(:, :, i), [ padsize, padsize]);
        eigfun = imresize(eigenfunction(:, :, i, xindex(1,j)), scalefactor);
        aerial=objectspectrum.*eigfun;%eigenfunction(:, :, i, j);
        
        FTaerial=fftshift(fft2(ifftshift(aerial)));
        I=I+x(j)*eigenvalue(i, xindex(1,j)).^2*abs(FTaerial).^2;
    end
end


if isreal(object)
    I_obj = normminmax(abs(object));
else
    I_obj = normminmax(angle(1-object));
end





switch qualitymethod
    case 'NRMSE'
        I_adj = normminmax(I);
        y = abs(NRMSE( I, I_obj ));
    case 'CPP'
        y = abs(1 - contrastperpixel( I ));
    case 'IMMSE'
        I_adj = normminmax(I);
        y = abs(immse( I, I_obj ));
    case 'SC'
        I_adj = normminmax(I);
        y = StructuralContent( I, I_obj );
    case 'PSNR'
        I_adj = normminmax(I);
        y = PeakSignaltoNoiseRatio( I, I_obj );
    case 'NCC'
        I_adj = normminmax(I);
        y = NormalizedCrossCorrelation( I, I_obj );
    case 'NAE'
        I_adj = normminmax(I);
        y = -NormalizedAbsoluteError( I, I_obj );
    case 'MSE'
        I_adj = normminmax(I);
        y = MeanSquareError( I, I_obj );
    case 'MD'
        I_adj = normminmax(I);
        y = MaximumDifference( I, I_obj );
    case 'AD'
        I_adj = normminmax(I);
        y = -AverageDifference( I, I_obj );
    case 'AMM'
        y = -AbsoluteMinMax(I);
    case 'STD'
        y = -std2(I);
    case 'CCP'
        I_adj = normminmax(I);
        y = -(max(max(I_adj))-min(min(I_adj)));
    otherwise
        
        error('wrong cost-function')
        
end








end



