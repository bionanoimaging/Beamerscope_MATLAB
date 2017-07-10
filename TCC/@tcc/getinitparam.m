function [ x0, illpattern ] = getinitparam( self, IOparams, INITparams, raw_object)
% function getinitparam = This function allows to get a first guess for the
% values of segments, archels or zernike coefficients resulting from the
% major diffraction orders of the objects FFT; 

% Possible patterning methods are 'Zernike', 'Archels', 'Segments'


N_mask = size(raw_object);

% Pad object for information which lays outside the usual frequency
% support; This is necessary for placing the pupils at positions outside
% the normal frequency grid
object_pad = raw_object;
object_pad = padarray(object_pad, [ size(object_pad, 1) size(object_pad, 2) ]);


% pad object with object copies
object = object_pad;
object = object + circshift(object, [0,N_mask]);
object = object + circshift(object, [N_mask,0]);
object = object + circshift(object, [N_mask,N_mask]);
object = object + circshift(object, [0,-N_mask]);
object = object + circshift(object, [-N_mask,0]);
object = object + circshift(object, [-N_mask,-N_mask]);
object = object + circshift(object, [-N_mask,N_mask]);
object = object + circshift(object, [N_mask,-N_mask]);

% Establish normalized coordinates.
%-----------------------------------------
v=self.x*(IOparams.NAo/IOparams.wavelength);
[self.vxx, self.vyy]=meshgrid(v);

% Prepare the normalized spatial-frequency grid.
% ---------------------------------------------------------
L=length(v);
vs=v(2)-v(1);                           % delta object coordinates
mc=1/(2*vs);                            % Spatial-frequency cutoff following nyquist .
self.m=linspace(-mc,mc,L);
self.m=self.m-min(abs(self.m));         % Ensure that DC position set to zero irrespective of numerical errors.
[self.mm, self.nn]=meshgrid(self.m);



% Obtain azimuthal and radial co-ordinates in the pupil plane.
%---------------------
[mtheta, mr]=cart2pol(self.mm,self.nn);
P2=(mr<=2); % Set all pixes outside of unit circle to zero for computation of Zernike polynomials.
P2_pad = padarray(P2, [ size(P2, 1) size(P2, 2) ],0);

% get zero filter
zerofilter_radius = 30/L;
P_zero = single(1-(mr<=zerofilter_radius));
P_zero_pad = padarray(P_zero, [ size(P_zero, 1) size(P_zero, 2) ],1);

% Initialize pupil of the system
%--------------------------------
Po=single(mr<=1);
Po_pad = padarray(Po, [ size(Po, 1) size(Po, 2) ],0);

% extend the grid to the size of the padded-object
% Po = padarray(Po, [ size(Po, 1) size(Po, 2) ]);

% get ojbect's spectrum (square of absulte value in logscale)
fft_spectr = sqrt(abs(fftshift(fft2(object)))+.001);
centre_fft = round(size(fft_spectr,2)/2);

% normalize and filter spectrum
fft_filtered = fft_spectr - min(min(fft_spectr));
fft_filtered = fft_filtered./max(max(fft_filtered));
fft_filtered = fft_filtered .*P2_pad.*padarray(P2, [ size(P2, 1) size(P2, 2) ],0); % cut out zero order; Pad pupil to fit extended frequency grid 


imagesc(sqrt(fft_filtered))
axis square
colormap gray


% %%%%%%idea: get order of descending diffraction-maxima in advance
% % get offset value
% zerothreshold= mean(mean(fft_filtered))*1;
% zerothreshold = 0.9;

% M=size(fft_filtered,2);
% N=numel(find(fft_filtered>zerothreshold)); % number of point sources
% [row,col]=find(fft_filtered>zerothreshold); % position of point sources
% illoverlap=zeros(size(Po,1), size(Po,2), N);
% minmaxorder = find(fft_filtered>zerothreshold);
% minmaxorder = [ minmaxorder row col];
% minmaxorder = sortrows(minmaxorder, -1);
%
% CP=zeros(size(object));
% for i=1:N%round(linspace(1,N,10),0)
%     CP=circshift(Po,round([minmaxorder(i,3)-(M+1)/2,minmaxorder(i,2)-(M+1)/2]));
%     illoverlap(:,:,i) = CP;
% end
% %toc


fft_spectr_temp = fft_filtered;
threshold_BO = mean(mean(fft_filtered));    % states the thrshold to differtiate between binary min/max diff order; for couting purposes

% initiliaze variables
iter=1;
clear mask_BO
n_difforders = 100;

while(n_difforders > 0 & iter < INITparams.maxIter)
    
    % find first maximum of the spectrum => max diffraction order
    [M,I] = max(fft_spectr_temp(:));
    [I_row, I_col] = ind2sub(size(fft_spectr_temp),I);
    
    % save coordinates for later use
    diffractionlist(1,iter) = I_row;
    diffractionlist(2,iter) = I_col;
    
    % count all the diffraction orders left; simply threshold the spectrum
    % and summ the values equals a rough guess of the number of diffraction
    % orders
    fft_spectr_thres = 1-im2bw(fft_spectr_temp, threshold_BO);
    n_difforders = sum(sum(fft_spectr_thres));
    
    % delete region with highest amount of BOs, surrounding results from
    % noise and has to get filtered anyway
    P_zero_temp = circshift(P_zero_pad, -[round(centre_fft - I_row) round(centre_fft - I_col)]);
    fft_spectr_temp = fft_spectr_temp.*P_zero_temp;
    
    % save current image of diffraction order for later superposition
    mask_BO(:,:,iter) = circshift(Po, -[round(centre_fft - I_row) round(centre_fft - I_col)]);
    
    
    
    %     imagesc(fft_spectr_temp)
    %     axis square
    %     drawnow
    %
    %     pause(0.1)
    
    iter = iter + 1;
end

imagesc(sum(mask_BO,3))
axis square

% sum all diffraction orders, normalize and filter according to max
% aperture of objective lens

overlap = sum(mask_BO,3).*P2;%.*Po;
overlap = (overlap/max(max(overlap)));

if(strcmp(IOparams.sourceshape,'Archels'))
    % TESTING Archels Method
    nMask = size(Po_pad,1);
    [illpattern] = diffr2archels(Po_pad, nMask, diffractionlist);
    
    illpattern = illpattern(1+N_mask(1):end-N_mask(2),1+N_mask(1):end-N_mask(2),:);
    
    for ( ii = 1:size( illpattern,3))
        x0(ii) = max(max(overlap.*illpattern(:,:,ii)));
    end
    
    x0 = double(x0);
    
%     subplot(121)
%     imagesc(sum(bsxfun(@times,illpattern,reshape(x0,[1 1 size(x0,2)])),3))
%     axis square
%     
%     subplot(122)
%     imagesc(overlap)
%     axis square
    
elseif(strcmp(IOparams.sourceshape,'MultiSegment'))
    % TESTING Archels Method
    
    
    r_segment = 3;      % number of radial segments (rings)
    n_segment = 12;     % number of circular segments
    
    IOparams.Zcoefficient = ones(1, r_segment*n_segment);
    
    for(i=0:round(r_segment*n_segment-1))
        
        % i-th segment in one of the annuli
        i_segment = rem(i, n_segment);
        
        if round(i/n_segment) == 0;
            NAc_i = 0;
            NAc_o = 0.33;
        elseif floor(i/n_segment) == 1;
            NAc_i = 0.33;
            NAc_o = 0.66;
        elseif floor(i/n_segment) == 2;
            NAc_i = 0.66;
            NAc_o = 1;
        end
        NAc_i = NAc_i * IOparams.NAc;
        NAc_o = NAc_o * IOparams.NAc;
        
        % Creating the annullar shape 0,1,2
        ith_seg=single(mr>=NAc_i/IOparams.NAo ... Inner radius.
            & mr<=NAc_o/IOparams.NAo ... Outer radius.
            );
        % scale rotational symmetric ramp 0..1
        segment_area = ((mtheta)/max(max(mtheta)) * round(n_segment/2)) + round(n_segment/2);
        
        % each segment results from the threshold of the grayvalues
        % filtered by the annular shape of the illumination sector
        % 0,1,2
        
        % this is due to aliasing of the pixelated source, otherwise
        % there will be a gap in the source shape
        if(i_segment == n_segment-1)
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1)*1.00001);
        else
            segment_area = double( segment_area >= i_segment & segment_area < (i_segment+1));
        end
        
        
        % get i-th segment and sum it up; weigh it with coefficient
        segment_area = segment_area.*ith_seg;
%         segment_area = padarray(segment_area, [ size(segment_area, 1) size(segment_area, 2) ],0);
        illpattern(:,:,i+1) = segment_area.*IOparams.Zcoefficient(1,i+1);
        
    end
    Ic = sum(illpattern, 3);
    
    % "Fit" diffraction orders to segments; give them their specific
    % weight!
    for ( ii = 1:size( illpattern,3) )
        x0(ii) = max(max(overlap.*illpattern(:,:,ii)));
    end
    
    x0 = double(x0);
    
elseif(strcmp(IOparams.sourceshape, 'Zernike'))
    
    source_raw = sum(mask_BO, 3);
    source_raw = source_raw./max(max(source_raw));
    source_raw = sqrt(source_raw);
    
    % h = fspecial('gaussian',11, 11);
    % I_res = imfilter(source, h, 'replicate');
    
    % fit zernike-polynom, get coefficients
    coeff = ZernikeCalc([1:IOparams.n_coefficients], source_raw);
    
    % draw result
    %ZernikeCalc([1:max_fit_order], coeff);
    % save result
    [result_sum illpattern] = ZernikeCalc([1:IOparams.n_coefficients], coeff, Po);
    centre_zern = round(size(result_sum,1)/2);
    
    source = result_sum( centre_zern-floor(size(Po, 1)/2):centre_zern+floor(size(Po, 1)/2), centre_zern-floor(size(Po, 1)/2):centre_zern+floor(size(Po, 1)/2));
    
    %source = source_raw( centre_fft-floor(size(raw_object, 1)/2):centre_fft+floor(size(raw_object, 1)/2), centre_fft-floor(size(raw_object, 1)/2):centre_fft+floor(size(raw_object, 1)/2));
    
    x0 = double(coeff');
    
else
    disp('Error! Wrong Source Shape!')
    return
end



end

