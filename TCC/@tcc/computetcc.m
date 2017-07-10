function tccResult=computetcc(self,sys,specimen,computeDevice)
% COMUTETCC: Computes images from the imaging system's transmission cross
% coefficient following Hopkins method
% functions contained in the microlith object and the specimen properties.
%   tcc=computetcc(MLOBJ, specimen,computeDevice) computes the TCC from the
%   transmission of the specimen/mask (in case of partially coherent
%   microscopy and lithography systems) or the density of fluorophores
%   (in case of incoherent microscopy) specified in specimen. The computed
%   image is returned and also saved within MLOBJ.
%
%   computeDevice = 'CPU', 'multiCPU', or 'GPU'.
%   'multiCPU' option computes images across the focus in parallel if the
%   MATLAB parallel computing toolbox is installed and a matlabpool is
%   available.
%
%
%   Written by Shalin Mehta, www.mshalin.com
%   License: GPL v3 or later.
%   Modified by Benedict Diederich

% This file is part of microlith package.
%
%     microlith is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     microlith is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with microlith.  If not, see <http://www.gnu.org/licenses/>.

if(~isa(self,'tcc'))
    error('This function requires an object of the microlith class.');
end

% Size and center of grid.
L=length(self.x);
DCAlongY=floor(1+0.5*size(self.Po,1));
DCAlongX=floor(1+0.5*size(self.Po,2));

disp('select function');
%integratePupil=@integratePupilPar; ATTENTION! NOT WORKING CORRECTLY
integratePupil=@integrateSOCSPAR;


% Read the pupils and allocate memory for the image.
%---------------------------------------------------
% We use single precision representation, because that runs the fastest on
% GPU. The round-off errors due to single-precision are not too expensive.

Tfun=single(self.Tfun);
tccResult=zeros(size(specimen,1),size(specimen,2),'single');


%% Compute image.
switch(self.config)
    
    case {'Brightfield','MultiSegment','Segment','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC','CustomPartiallyCoherent','DPC','Oblique','Dodt','DispersionStaining', 'Zernike'}
        
        %Find non-zero pixels on condenser and the amount of spectral shift
        %induced by them.
        %------------------------------------------
        Ic=single(self.Ic);
        [nshift,mshift, srcInt]=find(Ic);
        nshift=nshift-DCAlongY;
        mshift=mshift-DCAlongX;
        
        % Use single precision floating point to improve speed/memory.
        nshift=single(nshift);
        mshift=single(mshift);
        Ns=single(numel(srcInt));
        Ic=single(Ic);
        S=self.S;
               for idx=1:length(self.u)
                    %img(:,:,idx)=sumOverSource(objspec,Tfun(:,:,idx),Ic);
                    tccResult=integratePupil(Tfun,Ic,S);
                    
                end
        
end

% Apply the radiometric factor and assign to the object.
%-------------------------------------------------

% Radiometric factor that accounts for dependence of the image intensity
% on imaging NA and size of the sensor pixel. The effect of the condenser
% aperture relative to the imaging aperture is accounted for in the
% sumOverSource function.
% Our normalization is such that image of a 100nm^2 hole with an imaging
% and illumination NA of 1 is unity. For any other illumination NA, the
% intensity is proportional to ratio of the area of the condenser pupil to
% area of the objective pupil.


% ifftscale for grid-independent calculation of image in coherent
% and incoherent systems.
Np=numel(find(abs(self.Po)>1E-12));
ifftscale=L^2/Np;
%IFFT2 algorithm divides the input by numel(Po). The
%value at zero index of the output of IFFT2 is equal to the number of nonzero
%elements,i.e., numel(find(Po)). The above scale compensates
%for both and ensures that an image of a point produced by a clear
%circular aprture has a peak value of 1. This normalization allows us to
%compare images (apple to apple) computed over different grid sizes.

% Note that the numel(find(Po)) should be substituted by
% numel(find(Po*ObjectSpectrum)) when the object spectrum is
% finite and when the objectspectrum is shifted during partially coherent computation.
% However, for realistic specimens, the object spectrum almost
% always exceeds twice the support of Po and therefore the
% product of the pupil and shifted spectrum has the same support
% as the pupil.


PixSize=self.x(2)-self.x(1);


switch(self.config)
    case 'Coherent'
        RadiometricFactor=ifftscale*(PixSize/0.1)^2*(self.params.NAo)^2;
        tccResult=squeeze(RadiometricFactor*tccResult);
        
    case {'Brightfield','Segment','MultiSegment','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC','DPC', 'Oblique', 'Dodt', 'DispersionStaining', 'Zernike'}
        Sfactor=1/Np; % sum over soure 'brightens' the image by number of source points.
        % Normalization by the number of transparent points in imaging
        % pupil ensures that the intensity is proportional to Ns/Np.
        RadiometricFactor=ifftscale^2*Sfactor*(PixSize/0.1)^2*(self.params.NAo)^2;
        % For partially coherent imaging, ifft is applied on amplitude
        % image and then the result is squared. The radiometric factor
        % needs to take into account this squaring.
        tccResult=squeeze(RadiometricFactor*tccResult);
        
        
    case {'Fluorescence','Incoherent','Confocal'}
        RadiometricFactor=ifftscale*(PixSize/0.1)^2*(self.params.NAo)^2;
        tccResult=squeeze(real(RadiometricFactor*tccResult));
        % Fluorescence image has to be real and therefore objspec has to be conjugate-symmetric.
        % But, I don't use 'symmetric option of ifft2 because of the
        % preceding ifftshift.
end

self.img=tccResult;

%
%mkdir(strcat('result/',self.config,'/'))
%imwrite(imadjust((abs(specimen))), strcat('result/',self.config,'/',self.config, '_abs_NAo_', num2str(self.params.NAo), '-annulus_',num2str(self.params.annulus(1)),'_',num2str(self.params.annulus(2)), '.jpg'));
%imwrite(imadjust((angle(specimen))), strcat('result/',self.config,'/',self.config, '_angle_NAo_', num2str(self.params.NAo), '-annulus_',num2str(self.params.annulus(1)),'_',num2str(self.params.annulus(2)), '.jpg'));
%imwrite(imadjust(abs(img)), strcat('result/',self.config,'/',self.config,'_NAo_', num2str(self.params.NAo), '-annulus_',num2str(self.params.annulus(1)),'_',num2str(self.params.annulus(2)), '.jpg'));


end
