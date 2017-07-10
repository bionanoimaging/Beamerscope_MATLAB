%% decompose TCC into Kernels
function [h_fre h_vector h] =  decomposeTCC_SVD(TCC, order)
% this function decomposes a 2D matrix into their eigenvalues and eigenvectors 
%
%
% ----INPUT----
% TCC - is a 2D matrix calcullated from the overlap of the pupils
% order - specifies the number of kernels loading from the TCC
%
%
% ---OUTPUT----% 
% h_fre - 
% h_vector - 
%


N = sqrt(size(TCC,1));

disp('Calculating the SVD. Please wait...');
%%%%%%Singular value decomposition of the partially coherent imaging system%%%%%%
[U,S,V]=svd(TCC);   %Singular value decomposition

for(i=1:N)
    % frequency response kernel
    h_fre(:,:,i)=reshape(U(1:N^2,i:i),N,N);
    
    % spatial response kernel
    h(:,:,i)=(fftshift(ifft2(ifftshift(( h_fre(:,:,i))))));   %The impulse response of the first order approximation
end


sum_eigenvalue=(sum(sum(S)));   %The summation of the eigenvalues


%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%

% Get h-vector
h_vector = zeros(N, N);

disp('Calculating the h_vector. Please wait...');
for n_vect = 1:N
    for ii=1:N
        for jj=1:N
            h_vector((ii-1)*N+jj, n_vect)=h(ii,jj,n_vect);
        end
    end
end


% Get g-vector

disp('Calculating the g_1. Please wait...');
for n_vect=1:N
    for ii=1:N
        for jj=1:N
            g(ii,jj,n_vect)=h_vector((N-ii)*N+(N+1-jj),n_vect); %inverse vector
        end
    end
end


end

