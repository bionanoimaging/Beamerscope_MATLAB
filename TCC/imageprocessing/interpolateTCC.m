load('DPCNAc1_NAo1_rotAng180_xSize51.mat')

% testing image simulation
amplitude = im2double(imread('H:\Aufnahmen\AIDPC_ri_0_ro_392.5_1324_90_Method2_bud\resultnBright.tif'));
phase = im2double(imread('H:\Aufnahmen\AIDPC_ri_0_ro_392.5_1324_90_Method2_bud\resultDPCWaller_NAill0.95.tif'));
object_complex = amplitude .* exp(1i.*2*pi*(phase./max(max(phase))));


%%%%%%%%%%%%%%%%%%%%%%% temporary
DPCTCC = imrotate(DPCTCC, 0);

tic
% 2D->4D
TCC_x = reshapeTCC(DPCTCC);
t1 = toc

%Create a grid of points in  $R^4$. Then, pass the points through the function to create the sample values, V.
[x,y,z,t] = ndgrid(1:size(TCC_x,1),1:size(TCC_x,1),1:size(TCC_x,1),1:size(TCC_x,1));
t2 = toc

scaling_factor = .5:.1:1;
for(i=1:size(scaling_factor,2))
tic
%Now, create the query grid.
[xq,yq,zq,tq] = ndgrid(1:scaling_factor(i):size(TCC_x,1), 1:scaling_factor(i):size(TCC_x,1), 1:scaling_factor(i):size(TCC_x,1), 1:scaling_factor(i):size(TCC_x,1));
t3(i) = toc

tic
%Interpolate V at the query points.
TCC_interp = interpn(x,y,z,t,TCC_x,xq,yq,zq,tq,'cubic');
t4(i) = toc

tic
%reshape 4D->2D
TCC_interp_reshape = reshapeTCC4d2d( TCC_interp );
t5(i) = toc

object_complex_roi = object_complex(1:size(TCC_interp,1), 1:size(TCC_interp,1));

object_complex_roi = imresize(DPCsys.Ic,1/scaling_factor(i));
object_complex_roi = object_complex_roi(1:size(TCC_interp,1), 1:size(TCC_interp,1));
object_complex_roi = exp(1i*object_complex_roi*0.1/DPCparams.wavelength);

tic
disp('Calculating the aerial image. Please wait...')
%TCC_interp_reshape = DPCTCC;
aerial_freq = DPCsys.computeimagepar(TCC_interp_reshape, object_complex_roi, DPCparams, 'CPU');
t6(i) = toc
size_object(i) = size(object_complex_roi,1);


aerial=abs((ifft2((aerial_freq))))/((size(aerial_freq,2))^2);

figure
subplot(1,3,1)
(imagesc(aerial))
axis square
colormap gray

subplot(1,3,2)
(imagesc(abs(object_complex_roi)))
axis square
colormap gray

subplot(1,3,3)
(imagesc(angle(object_complex_roi)))
axis square
colormap gray

end
% % Visualization
% for(i=1:size(TCC_interp,1))
% figure(1)
% imagesc(squeeze((TCC_interp(round(size(TCC_interp,1)/2),i,:,:))))
% axis square
% colormap jet
% drawnow
% pause(0.1);
% end

time = [t3; t4; t5; t6]

figure
bar(size_object, time')
xlabel size_object
ylabel 'time in s'
xlabel 'objectsize in px'
legend('Calculate Upsampled Grid', 'Interpolate TCC', 'Reshape TCC 4D-2D', 'Calculate Image', 'Location','northwest')

figure
plot(size_object, t6)
hold on
plot(size_object, t5)
plot(size_object, t4)
plot(size_object, t3)
ylabel 'time in s'
xlabel 'objectsize in px'
hold off

save('./Data/TimeInterpolationExperiment.mat', 'time', 'size_object')


matlab2tikz('./tikz/bar_chart_computationaltime.tex', 'height', '\fheight', 'width', '\fwidth' )


% %Create a movie to show the results.
% figure('renderer','zbuffer');
% nframes = size(tq, 4);
% for j = 1:nframes
%    slice(yq(:,:,:,j),xq(:,:,:,j),zq(:,:,:,j),...
%          TCC_x(:,:,:,j),0,0,0);
%    caxis auto;
%    colormap jet
%    pause(0.1)
%    drawnow
% end
