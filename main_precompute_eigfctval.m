%% Generate a dataset with complex objects and its corresponding optimized
% illumination shapes using the TCC
% 
% The software is for private use only and gives no guarrantee, that it's
% working as it should! 
% 
%
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.


%% ANDROID-SECTION
% clear all
close all


%% make sure all components are loaded correctly and path is correct

% WINDOWS
root_path = 'C:\Users\Bene\Dropbox\Dokumente\OpticalDesign\Beamerscope\MATLAB\';

% MAC
root_path = '/Users/Bene/Dropbox/Dokumente/Promotion/PROJECTS/Beamerscope/MATLAB/';

% add path of all libraries used in this process (optimization Toolbox from
% MATLAB is required!
addpath(genpath(strcat(root_path, 'Beamerscope_MATLAB')))
addpath(genpath(root_path))
cd([root_path, 'Beamerscope_MATLAB'])




%% System parameters (here the ones from the cellphone microscope)
% Set parameters of the acquisition process (in mum)
wavelength=0.530;       % systems center-wavelength
NAo = 0.25;             % aperture of microscopes objective
NAc = 0.5;              % maximum condenser aperture

pixel_eff = 0.6;        % representing effective pixelsize of the camera
pixel_n = 128;           % number of simulated pixel

n_immers = 1.00;        % refractive index immmersion medium (air)
n_obj = 1.20;           % refractive index object

%% Set parameters of the optimization process
sourceshape = 'MultiSegment';    % base for source shape optimization ('Zernike', 'MultiSegment', 'Archels')
n_coefficients = 48;             % number of segments/zernikes used for optimizaton



%%%%%%%%Programm starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate spatial/frequency grid

simRadius=((pixel_n-1)*pixel_eff)/2;      % radius of simulated region in microns
xsim=-simRadius:pixel_eff:simRadius;      % create grid for simulation
zsim= 0;                                  % simulated z-position (not yet implemented)
usim=0;                                   %

kernel_limit = 50;                        % limit the number of eigenfunctions for simulation (default = 50)


%% Create paramter's object; IO = Illumination Optimization
IOparams.wavelength=wavelength;
IOparams.NAo=NAo;
IOparams.NAc=NAc;
IOparams.pixel_eff = pixel_eff;
IOparams.pixel_n = pixel_n;
IOparams.path = path;
IOparams.sourceshape  = sourceshape;
IOparams.n_coefficients = n_coefficients;
IOparams.kernel_limit = kernel_limit;
IOparams.illpattern = 0;
IOparams.xsim = xsim;

%% precalculate test-TCC and eigenvalues for desired parameterset (according to light-source shape etc.)
[ IOsys ] = precalculatetcc( IOparams, xsim, zsim );

%% assign eigenvalues/eigenfunctions
eigenvalue = IOsys.eigenvalue;
eigenfunction = IOsys.eigenfunction;

%% fftshift the kernels before saving - Tensorflow requiers that! 
eigenfunction_shift = {};
index = 1;
for j = 1:size(eigenfunction, 4)
    for i= 1:size(eigenfunction, 3)
        eigenfunction_shift{index}=ifftshift(eigenfunction(:,:,i,j));
        index = index + 1;
    end
end

eigenfunction_shift = cat(4, eigenfunction_shift{:});

%%
save('/Users/Bene/Dropbox/Dokumente/Promotion/PROJECTS/Beamerscope/Python/Beamerscope_TENSORFLOW/Beamerscope_IllOpt/tf_illoptdata.mat', 'eigenvalue', 'eigenfunction', 'ill_method', '-v7.3')




