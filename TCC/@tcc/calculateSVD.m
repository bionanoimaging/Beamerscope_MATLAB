function [eigenfunction eigenvalue_return] = calculateSVD(self)
% This function decompose the TCC into its eigenvectors and eigenfunctions
% according to the paper:
% Aerial Image Simulation for partial coherent system with programming
% development in MATLAB
%
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.
%
%
%
% input:
% - TCC - 2D TCC representation; each TCC is stacked into one row
% - M, N - Dimensions of the TCC
%
%
% output:
% - eigenfunction
% - eigenvalue

if(~isa(self,'tcc'))
    error('This function requires an object of the microlith class.');
end

M = self.M;
N = self.N;

%disp('%%%%%%Calculate SVD and Eigenfunctions%%%%%%%')
%tic
% Apply SVD to singular matrix P we will get eigenvalues and eigenfunctions:
try
    [U, eigenvalue, V]=svd(self.TCC,'econ'); % Singular value decomposition
catch
    warning('Problem using SVD function.');
    
    self.TCC(isnan(self.TCC))=0;
    self.TCC(isinf(self.TCC))=0;

    [U, eigenvalue, V]=svd(self.TCC,'econ'); % Singular value decomposition
end


% To get those eigenfunctions:
eigenfunction = zeros(M,M,N);
self.eigenvalue = zeros(1,N);
for ii=1:N
    if(ii<=size(V,2))
        eigenfunction(:, :, ii)=reshape(V(:, ii),M,M)'; % plot eigenfunction
        
        %  imagesc(eigenfunction(:, :, i))
        %  axis square
        %  colormap
        %  drawnow
        self.eigenvalue(ii) = eigenvalue(ii,ii);
    end
end
eigenvalue_return = self.eigenvalue';
self.eigenfunction = eigenfunction;


%toc