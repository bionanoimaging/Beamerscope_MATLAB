function computesys(self,config,params)
% COMPUTESYS: Computes the transfer function of the imaging path and
% the intensity distribution of the illumination aperture for given
% microscopy/lithography system.
%   COMPUTESYS(MLOBJ,config,params) calculates the transfer properties
%   using the config string and parameters specified in the params
%   structure.
%
%   Valid string-values for config are:
%   'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC',
%   'Fluorescence','Incoherent','Confocal','Coherent', 'Custom'...
%   Following table summarizes the parameters that can be present in the params structure
%   and for which configurations they are used.
%     NAo           Numerical aperture of the imaging path (all configurations)
%     wavelength    Wavelength (in um) that affects the optical resolution
%     (all configurations). If wavelength is a vector, polychromatic
%     illumination is assumed.
%     NAc           Numerical aperture of the condenser (only for transmitted light methods).
%     annulus       For dark-field and phase-contrast systems the annulus
%                   parameter has to be a vector [InnerNA OuterNA].
%     phasering     For phase-contrast system phasering parameter has to be a
%                   vector [InnerNA OuterNA Absorption PhaseDelay]
%     annulus_kond  amount how much the aperture in the focal plane of the condensor is open
%     annlusu_konj  amount how much the aperture in the focal plane of the microscope lens is open
%     shear         For DIC,DIC-Preza, and PlasDIC systems,
%                   shear specified in units of wavelength/NAo.
%     bias          For DIC,DIC-Preza, and PlasDIC systems,
%                   the phase-bias specified in degrees such that bias of
%                   90 leads to brightfield contrast.
%     shearangle    For DIC,DIC-Preza, and PlasDIC systems,
%                   the direction of shear specified in degrees.
%     CameraSensitivity Sensitivity of the camera as a function of
%                       wavelength.
%     SourceSpectrum    Spectral distribuition of the source as a function
%                       of wavelength.
%
%     astigmatism   Astigmatism specified as coefficients of Zernike polynomials (n,m)=(2,-2) and
%     (2,2)
%     coma          Coma specified as coefficeients of Zernike polynomials
%     (n,m)= (3,1) and (3,-1).
%     spherical     Spherical aberration specified as coefficent for Zernike polynomial (n,m)=(4,0).
%     ChromaticShiftX
%     ChromaticShiftY chromatic shift across wavelengths mentioned in wavelength.
%     ChromaticScale chromatic scaling (beyond the scaling caused by diffraction)

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


%   The optical units used in our simulation are units in Born & Wolf  divided by 2pi. In
%   the units used by Born & Wolf, the first zero of airy disk occurs at
%   0.61*2pi=3.83 and the first on-axis minima occurs at normalized defocus of 4pi.
%   In the optical units used for simulation, the first zero of airy disk occurs at 0.61 and the
%   axial minima at 2.
%

% Establish normalized coordinates.
%-----------------------------------------

v=self.x*(params.NAo/params.wavelength);
[self.vxx, self.vyy]=meshgrid(v);

% Prepare the normalized spatial-frequency grid.
% ---------------------------------------------------------
L=length(v);
vs=v(2)-v(1); % delta object coordinates
mc=1/(2*vs); % Spatial-frequency cutoff following nyquist .
self.m=linspace(-mc,mc,L);
self.m=self.m-min(abs(self.m)); % Ensure that DC position set to zero irrespective of numerical errors.
[self.mm, self.nn]=meshgrid(self.m);


% Obtain azimuthal and radial co-ordinates in the pupil plane.
%---------------------
[mtheta, mr]=cart2pol(self.mm,self.nn);
mrUnitCircle=mr.*(mr<=1); % Set all pixes outside of unit circle to zero for computation of Zernike polynomials.

% Initialize pupils of the system
%--------------------------------

Po=single(mr<=1);



% Now set-up the pupils according to configuration
%-----------------------

Tfun=complex(zeros(size(Po),'single'),zeros(size(Po),'single'));

%Tfun=self.Tfun;
switch(config)
    case {'Brightfield'}
        S=params.NAc/params.NAo
        Ic=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.
        Tfun = Po;
        
        
    case {'Archels'}
        S=params.NAc/params.NAo;
        Ic = zeros(size(Po));
        
        % get i-th segment and sum it up; weigh it with coefficient
        Ic = sum(bsxfun(@times,params.illpattern,reshape(params.Zcoefficient,[1 1 size(params.Zcoefficient,2)])),3);
        
        Tfun = Po;
        
        
    
        
    case {'MultiSegment'}
        S=params.NAc/params.NAo;
        n_segment = 12;      % number of colour segments
        n_rings = 4;
        Ic = zeros(size(Po));
        
        for(i=0:size(params.Zcoefficient,2)-1)
            
            Isegment = zeros(size(Po));
            
            % i-th segment in one of the annuli
            i_segment = rem(i, n_segment);
            
            if floor(i/n_segment) == 0;
                NAc_i = 0;
                NAc_o = 0.25;
            elseif floor(i/n_segment) == 1;
                NAc_i = 0.25;
                NAc_o = 0.5;
            elseif floor(i/n_segment) == 2;
                NAc_i = 0.5;
                NAc_o = .75;
            elseif floor(i/n_segment) == 3;
                NAc_i = 0.75;
                NAc_o = 1;
            end
            NAc_i = NAc_i * params.NAc;
            NAc_o = NAc_o * params.NAc;
            
            % Creating the annullar shape 0,1,2,3
            Isegment=single(mr>=NAc_i/params.NAo ... Inner radius.
                & mr<=NAc_o/params.NAo ... Outer radius.
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
            segment_area = segment_area.*Isegment;
            Isegment = segment_area.*params.Zcoefficient(1,i+1);
            
            Ic = Ic + Isegment;
            
            
        end
        
        Tfun = Po;
        
        
    case 'Darkfield'
        % Set the condenser annulus.
        Ic=single(mr>=params.NAci/params.NAo ... Inner radius.
            & mr<=params.NAc/params.NAo ... Outer radius.
            );
        
        Tfun=Po;
        S=params.NAc/params.NAo;
        
        
    case 'PhaseContrast'
        % Set the condenser and objective annuli.
        % The pupil co-ordinate is normalized with respect to NA
        % of the objective.
        
        Ic=single(mr>=params.annulus(1)/params.NAo ... Inner radius.
            & mr<=params.annulus(2)/params.NAo ... Outer radius.
            );
        
        Ic=single(mr<=params.annulus(2)/params.NAo);
        Ic(size(Ic,1)/2:end, :) = 0;
        
        phasering=params.phasering(3)* exp(1i*params.phasering(4)) *...
            single(mr>=params.phasering(1)/params.NAo &...
            mr<=params.phasering(2)/params.NAo);

        
        phasering = (mr<=params.phasering(2)/params.NAo);
        phasemask = ones(size(phasering));
        phasemask(size(phasering,1)/2:end, :)=params.phasering(3)* exp(1i*params.phasering(4)); 
        
        phasering = phasering.*phasemask;
        Po = single(phasering);
        
        
        %[InnerNA OuterNA Absorption PhaseDelay]
        %Direct multiplication by phasering sets transmission
        %outside phasering to zero. Therefore, create logical
        %mask to place the phasering in the objective pupil.
        
        %Po = zeros(size(phasering));
        %phaseringmask=(phasering~=0);
        %Po(phaseringmask)=phasering(phaseringmask);
        
        %Po = zeros(size(Po));
        %Po(size(Po,1)/2:end, size(Po,2)/2:end) = 1;
        
        Tfun=Po;
        S=params.NAc/params.NAo;
        
        
    case 'DispersionStaining'
        if(~ ((isfield(params,'annulus') || isfield(params,'fieldstop'))) )
            error('For DispersionStaining system, the params structure must have the ''annulus'' or the ''fieldstop''');
        end
        
        S=params.NAc/params.NAo;
        
        if(isfield(params,'annulus') )
            % Set the condenser and objective annuli.
            % The pupil co-ordinate is normalized with respect to NA
            % of the objective.
            
            Ic=single(mr>=params.annulus(1)/params.NAo ... Inner radius.
                & mr<=params.annulus(2)/params.NAo ... Outer radius.
                );
            
        elseif(isfield(params, 'fieldstop'))
            
            %moPupil=single(mr>=params.fieldstop(1)/params.NAo ... Inner radius.
            % & mr<=params.fieldstop(2)/params.NAo ... Outer radius.
            %  );
            
            moPupil=1-single(mr<=params.fieldstop(2)/params.NAo);
            
            Ic=single(mr<=S);
            Po = Po.*moPupil;
        end
        
        
        parfor idx=1:length(u)
            Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
        end
        
    case {'Zernike'}
        S=params.NAc/params.NAo;
        Ic_mask=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.
        Tfun = Po;
        
        [Ic temp] = (ZernikeCalc(1:size(params.Zcoefficient,2), params.Zcoefficient', Ic_mask, [], 'standard'));
        Ic = Ic-min(min(Ic));
        Ic = Ic/max(max(Ic));
        Ic = Ic.*Ic_mask;
        
    case 'DPC'
        
        %         if(~ ( isfield(params,'NAc') && isfield(params,'NAo') && isfield(params,'rot_cond') && isfield(params,'annulus_cond')) )
        %             error('For phase-gradient system, the params structure must have the ''annulus_kond'' (1x1) and the ''annulus_konj'' (1x2) vectors.');
        %         end
        % Set the condenser and objective annuli.
        % The pupil co-ordinate is normalized with respect to NA
        % of the objective.
        
        
        S=params.NAc/params.NAo;
        aperture_condensor=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.
        
        
        radiusAperture = round((S-0.5)*size(aperture_condensor,2));
        centerIndex = round(size(aperture_condensor,2)/2);         %# Get the center index for the cols
        cutArea_condensor = uint8(round(centerIndex + .5*radiusAperture*(.5-params.annulus_cond)));
        
        aperture_condensor(:, cutArea_condensor:end) = 0;
        aperture_condensor = imrotate(aperture_condensor, params.rot_cond ,'bilinear','crop');
        
        Ic = aperture_condensor;
        Tfun=Po;
        
        % %        Test if changing MO pupil instead of PC equals in the same result
        %         Ic = Po;
        %         Tfun= aperture_condensor;
        %
        
    case 'Dodt'
        % Set the condenser annulus.
        if(~isfield(params,'annulus')|| ~isfield(params,'rot_angle')|| ~isfield(params,'gaussSize')|| ~isfield(params,'gaussSpread'))
            error('For Dodt system, the params structure must have the ''annulus'' (a 1x2 array with innerNA and outerNA) and ''rot_angle''');
        end
        
        Ic = zeros(size(mr));
        Temp = single(mr>=params.annulus(1)/params.NAo ... Inner radius.
            & mr<=params.annulus(2)/params.NAo ... Outer radius.
            );
        Ic(round(size(Ic,2)/2):end,round(size(Ic,2)/2):end)=Temp(round(size(Ic,2)/2):end,round(size(Ic,2)/2):end);
        
        Ic = imrotate(Ic, params.rot_angle ,'bilinear','crop');
        
        
        
        
        f_gauss = fspecial('gaussian', params.gaussSize, params.gaussSpread );
        Ic = filter2(f_gauss, Ic);
        
        
        parfor idx=1:length(u)
            Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
        end
        
        
    case 'Oblique'
        
        if(~ (isfield(params,'annulus_cond') || isfield(params,'rot_angle')) )
            error('For oblique system, the params structure must have the ''annulus_kond'' (1x1), the ''annulus_konj'' (1x2) and ''rot_angle'' vectors.');
        end
        % Set the condenser and objective annuli.
        % The pupil co-ordinate is normalized with respect to NA
        % of the objective.
        
        
        S=params.NAc/params.NAo;
        aperture_condensor=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.
        
        
        radiusAperture = round((S-0.5)*size(aperture_condensor,2));
        centerIndex = round(size(aperture_condensor,2)/2);         %# Get the center index for the cols
        cutArea_condensor = (round(centerIndex + .5*radiusAperture*(.5-params.annulus_cond)));
        
        aperture_condensor(:, cutArea_condensor:end) = 0;  %# Set the right half to the value 0 (of the same type as A)
        aperture_condensor = imrotate(aperture_condensor, params.rot_angle ,'bilinear','crop');
        
        
        aperture_conjugate =single(mr<=S);
        
        
        Ic = aperture_condensor;
        
        
        
        %Direct multiplication by phasering sets transmission
        %outside phasering to zero. Therefore, create logical
        %mask to place the phasering in the objective pupil.
        
        aperture_mask=(aperture_conjugate~=0);
        %Po(aperture_mask)=aperture_konjugate(aperture_mask);
        Po = Po .* aperture_conjugate;
        
        parfor idx=1:length(u)
            Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
        end
        
        
    case 'Zernike'
        
        S=params.NAc/params.NAo;
        Ic=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.
        Tfun = Po;
        
        % Po(Po<0)=1; %For pixels that are entirely within pupil set the transmittance to 1.
        % Po(Po>1)=0; %For pixels that are entirely outside pupil set the transmittance to 0.
        
        % Model aberrations (astigmastism, coma, spherical) using Zernike
        % polynomials
        %---------------------------------
        if(isfield(params,'astigmatism'))
            Pastigm=params.astigmatism(1)*mrUnitCircle.^2.*cos(2*mtheta) +...
                params.astigmatism(2)*mrUnitCircle.^2.*sin(2*mtheta);
        else
            Pastigm=zeros(size(Ic));
        end
        
        if(isfield(params,'coma'))
            Pcoma=params.coma(1)*(3*mrUnitCircle.^3-2*mrUnitCircle).*cos(mtheta)+...
                params.coma(1)*(3*mrUnitCircle.^3-2*mrUnitCircle).*sin(mtheta);
        else
            Pcoma=zeros(size(Ic));
        end
        
        if(isfield(params,'spherical'))
            Pspherical=params.spherical*(6*mrUnitCircle.^4-6*mrUnitCircle.^2+1);
        else
            Pspherical=zeros(size(Ic));
        end
        
        ICAberration=exp(1i*Pastigm).*exp(1i*Pcoma).*exp(1i*Pspherical);
        Ic=angle(Ic.*single(ICAberration));
        Ic=Ic-min(min(Ic));
        disp('Zernike')
        
end

self.config=config;
self.params=params;
self.Ic=Ic;
self.Po=Po;
self.Tfun=Tfun;
self.v=v;
self.S=S;   %coherence factor
end

function otf=ctf2otf(ctf)
% OTF=CTF2OTF(CTF) computes the optical transfer function (OTF) from
% the coherent transfer function (CTF) of the imaging system.
% OTF is the auto-correlation of CTF, which we compute using the FFTs.
% The method generates the OTF on the same grid as CTF and hence the
% OTF has twice the cut-off frequency as that of CTF.
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
%IFFT2 algorithm divides the input by numel(Po). The
%value at zero index of the output of IFFT2 is equal to the number of nonzero
%elements,i.e., numel(find(Po)). The above scale compensates
%for both and ensures that an image of a point produced by a clear
%circular aprture has a peak value of 1. This normalization allows us to
%compare images (apple to apple) computed over different grid sizes.

apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;

% The ifftscale should be applied to ipsf, because if applied to apsf; it
% is applied twice.

otf=fftshift(fft2(ifftshift(ipsf)));

end
