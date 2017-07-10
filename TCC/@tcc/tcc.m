classdef tcc <handle
    % TCC class implements calculation of the transmission cross coefficient of
    % a microscopy and optical lithography systems.
    % While the package allows simulation of general microscopy system,
    % following imaging methods can be simulated with simple interface.
    % * Brightfield
    % * Darkfield
    % * Phase-contrast
    % * Differential Interference Contrast
    % * Widefield fluorescence
    % * Confocal fluorescence
    % * Coherent imaging system
    %
    % This work is majorly derived from Shalin Methas "Microlith"-Toolbox and
    % extends its licensing rules!
    %
    %   The image computation is carried out in three steps as follows.
    %
    %   Step 1: Create a TCC object (allocates the required memory).
    %   USAGE:  TCCOBJ= TCC(x,z);
    %   Using a vector for z parameter simulates 3D image.
    %
    %   Step 2: Compute transfer functions that represent the imaging
    %   system. USAGE:  TCCOBJ.computesys(config,parameters) OR
    %           computesys(TCCOBJ,config,parameters)
    %
    %   Step 3: Compute the 2D/3D image using specimen's 2D transmission
    %   or fluorescence described over a square grid defined by x. USAGE:
    %   MLOBJ.computeimage(specimen) OR
    %           computeimage(MLOBJ,specimen)
    %
    %   The user can simulate a general microscopy or lithography system by
    %   substituting the objective and condenser pupils after the step-2 and
    %   before the step-3.
    %
    %   The mutual intensity in the focal region (only for coherent/partially
    %   coherent systems) can be computed with TCCOBJ.computecoherence() OR
    %   computecoherence(MLOBJ)
    %
    %
    %
    % A TCC object follows handle semantics, i.e., the functions operating
    % on the object do not modify the copy, but the present instance of the object.
    %
    %
    %   Written by Shalin Mehta, www.mshalin.com
    %   License: GPL v3 or later.
    %   Extended by Benedict Diederich
    
    % This file is part of tcc package.
    %
    %     tcc/TCC is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     tcc is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with tcc.  If not, see <http://www.gnu.org/licenses/>.
    
    properties
        
        Tfun        % The transfer function of the imaging path.
        % Transmission of the objective pupil
        % for non-fluorescent microscopy or lithography
        % configurations. OTF for fluorescence
        % configurations.
        Po          % Objective pupil.
        Ic          % Condenser intensity distribution.
        % Relevant only for transmission/reflection
        % microscopy and lithography configuration.
        S           % Coherence Factor S=NAc/NAo;
        
        x          % Physical support in the specimen plane (in um).
        v          % Support in the specimen plane (in  normalized units).
        m          % Support in the frequency domain (in  normalized units).
        mu         % Support in the frequency domain along the optical axis.
        vxx        % vxx and vyy define the spatial grid (in normalized units).
        vyy
        z           % Focal distances (z=o implies the focal plane) (in um).
        u           % Focal distances (in normalized units).
        mm         % mm and nn define the spatial-frequency grid (in normalized units).
        nn
        specimen    % 2D transmission/reflection for microscopy and lithography system.
        % 2D fluorophore density for fluorescence imaging
        % systems.
        img         % 3D image (with lateral support the same as x and the axial support same as z).
        Zcoefficient% Vector with Zernike Coefficients
        config % String that identifies the type of imaging system:
        %'BF' Brightfield, 'DF' Darkfield, 'PC' Phase-contrast, 'FL1P' Standard (1-photon) Fluorescence, 'FL2P' Two photon fluorescence, 'FLC' Confocal fluoresence.
        params %  Parameters of the system as well as aberrations of the imaging pupil specified as parameter value pairs.
        %             NAo         % Numerical aperture of the imaging path
        %             NAc         % Numerical aperture of the condenser (relevant only for transmitted light methods).
        %             nImm          % Refractive index of the objective immersion medium.
        %             nEmbb        % Refractive index of the specimen embedding medium.
        %             wavelength      % Wavelength that governs the resolution.
        
        M       % Size of TCC Matrix
        N       % Size of TCC Matrix
        TCC     % TCC Matrix
        eigenfunction   % eigenfunction of the decomposed TCC
        eigenvalue      % eigenvalues of the decomposed TCC
        
        indexSegment    % i-th segment in segmented illumination
        intensity       % intensity of segemtn
        
    end
    
    methods
        function self=tcc(x,z)
            
            %   MIC = MICROLITH(x,z) allocates memory for an
            %    object that represents a microscopy or lithography imaging
            %    system.
            %   Vector x defines support along both X and Y axis (i.e., we
            %   assume imaging over a square lateral grid), and focal distances
            %   are specified by vector z (z=0 implies the focal plane). x, z,
            %   and wavelength of imaging are assumed to be expressed in the
            %   same units (um works well for optical microscopy).
            %
            %   Written by Shalin Mehta, www.mshalin.com
            %   License: GPL v3 or later.
            
            
            
            
            % Allocate the spatial grid.
            % -------------------
            self.x=x;
            self.v=zeros(size(x),'double');
            self.z=z;
            self.u=zeros(size(z),'double');
            self.vxx=zeros(length(x),length(x),'double');
            self.vyy=zeros(length(x),length(x),'double');
            
            % Allocate the spatial frequency grid
            % ---------------------
            self.m=zeros(size(x),'double');
            self.mm=zeros(length(x),length(x),'double');
            self.nn=zeros(length(x),length(x),'double');
            
            % Allocate memory for the system pupils, specimen, and the
            % image.
            %-----------------------------------
            self.specimen=zeros(length(x),length(x),'double');
            self.img=zeros(length(x),length(x),length(z),'double');
            self.Tfun=zeros(length(x),length(x),'double');
            self.Ic=zeros(length(x),length(x),'double');
            
        end
        
        computesys(self,config,params); % Compute transfer function of the imaging path based on the parameters.
        tccResult=computetcc(self,params,computeDevice);
        img = computeimage(self, TCCmat, objspec, params, computeDevice);
        
        % Major Calculation is done here
        [TCC M N] = calculateTCC(self);
        
        % Decomposing the TCC into its eigenvectors and eigenvalues
        [eigenfunction eigenvalue] = calculateSVD(self);
        
        % Calculate the aerial image from the TCC and object
        [I] = calculateImage(self, object);
        
        % Optimize Illumination source
        [ x_opt,fval,output ] = optillsource( self, IOparams, object, x0 )
        
        % get inital parameters (archels/Siedentopf principle)
        [ x0, illpattern ] = getinitparam( self, IOparams, INITparams, raw_object)
    end
    
end
