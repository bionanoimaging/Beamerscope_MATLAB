%% decompose TCC into Kernels
% clear all;
% matlabpool(4)
addpath(genpath('../TCC/'))
%load('C:\Users\Bene\Dropbox\Dokumente\Masterarbeit\Matlab\TCC\Data\ZernikeNAc1_NAo1_ZernikePolynom3_xSize51.mat')
% load('..\TCC\Data\DPCNAc1_NAo1_rotAng0_xSize51.mat')

TCC4D = TCC_x;
TCC = imrotate(TCC_SOCS,90);%imrotate(reshapeTCC4d2d(TCC4D),0);

N_mask = size(TCC4D, 1);
N_mask = 31;
figure
imagesc(TCC)
axis square


tic
disp('Calculating the SVD. Please wait...');
%%%%%%Singular value decomposition of the partially coherent imaging system%%%%%%
[U,S,V]=svd(TCC);   %Singular value decomposition

clear h_fre h
for(i=1:N_mask)
    % frequency response kernel
    h_fre(:,:,i)=reshape(U(1:N_mask^2,i:i),N_mask,N_mask);
    
    % spatial response kernel
    h(:,:,i)=(fftshift(ifft2(ifftshift(( h_fre(:,:,i))))));   %The impulse response of the first order approximation
end
toc

sum_eigenvalue=(sum(sum(S)));   %The summation of the eigenvalues


%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%

% Get h-vector
h_vector = zeros(N_mask, N_mask);
tic
disp('Calculating the h_vector. Please wait...');
for n_vect = 1:N_mask
    for ii=1:N_mask
        for jj=1:N_mask
            h_vector((ii-1)*N_mask+jj, n_vect)=h(ii,jj,n_vect);
        end
    end
end
toc

% Get g-vector
tic
disp('Calculating the g_1. Please wait...');
for n_vect=1:N_mask
    for ii=1:N_mask
        for jj=1:N_mask
            g(ii,jj,n_vect)=h_vector((N_mask-ii)*N_mask+(N_mask+1-jj),n_vect); %inverse vector
        end
    end
end
toc

%%%%%The desired output pattern%%%%%%
pz=zeros(N_mask,N_mask);
for ii=16:35
    for jj=13:23
        pz(ii,jj)=1;
    end
end
for ii=16:35
    for jj=29:39
        pz(ii,jj)=1;
    end
end

pz = exp(1i*pz*0.003);
pz = object(11:41,11:41);
%%%%%%The initialization of \theta, where r=\theta%%%%%%
r=ones(N_mask,N_mask)*pi/2;
for ii=16:35
    for jj=13:23
        r(ii,jj)=pi/5;
    end
end
for ii=16:35
    for jj=29:39
        r(ii,jj)=pi/5;
    end
end



disp('Calculating the SVD-Image. Please wait...');
tic
%%%%%%PSM optimization in partially coherent imaging system%%%%%%
clear aerial_svd
aerial_svd_i = zeros(N_mask, N_mask);
for(number_kernel=1:2)
h_res = h(:,:,number_kernel);
aerial_svd(:,:,number_kernel) = (  (imfilter(double(pz),h_res)));
aerial_svd_i = aerial_svd_i+aerial_svd(:,:,number_kernel);
imagesc(abs(h_res))
axis square
colormap gray
pause(0.2)
end
toc
aerial_svd_i = abs(aerial_svd_i).^2;
kernel = [-1, -1, -1, -1, 8, -1, -1, -1]/8;
diffImage = conv2(double(aerial_svd_i), kernel, 'same');
cpp = mean2(diffImage)

figure
subplot(132)
imagesc(abs(aerial_svd_i.^2));
axis on;
axis square
colormap gray
title('Output of desired pattern - SVD-Plot');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Output pattern of desired pattern%%%%%%
disp(strcat('Start computing the image using classic TCC implementation. Please wait...'));

%% mathwork solution
tic

pz_fre=(fftshift(fft2(pz)));

[index_x, index_y] = ndgrid(1:N_mask^2);

[index_1,index_2] =  ind2sub([N_mask,N_mask],  index_x);
[index_3,index_4] =  ind2sub([N_mask,N_mask],  index_y);

index_m = (mod(index_1-index_3,N_mask) + 1).';
index_n = (mod(index_2-index_4,N_mask) + 1).';

vals=bsxfun(@times, TCC',pz_fre(:));
vals=bsxfun(@times,vals,pz_fre(:)');

idx=(vals~=0);

subs=[index_m(idx), index_n(idx)];
vals=vals(idx);

aerial_fre = accumarray(subs,vals,[N_mask,N_mask]);    
toc


aerial=abs(ifft2(aerial_fre))/((N_mask)^2)./N_mask;

kernel = [-1, -1, -1, -1, 8, -1, -1, -1]/8;
diffImage = conv2(double(aerial), kernel, 'same');
cpp = mean2(diffImage)
%% display aerial image

aerial_fre=computeimageparXX(TCC, pz);

subplot(131)
imagesc(aerial);
axis on;
axis square
colormap gray
title('Output of desired pattern');

aerial = aerial./max(max(aerial));
aerial_svd_i = aerial_svd_i./max(max(aerial_svd_i));

subplot(133)
plot(1:N_mask,aerial(round(N_mask/2),:));
hold on
plot(1:N_mask,aerial_svd_i(round(N_mask/2),:));
hold off
axis on;
axis square
colormap gray
title('Output of desired pattern');




% %% display aerial image
% subplot(133)
% imagesc(pz-aerial);
% axis on;
% axis square
% colormap gray
% title('Error of desired pattern');
