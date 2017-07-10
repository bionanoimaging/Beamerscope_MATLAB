%% ANDROID-SECTION
% clear all
close all


% make sure all components are loaded correctly and path is correct

%root_path = 'H:\'; % WINDOWS
root_path = '/Users/Bene/Dropbox/Dokumente/OpticalDesign/Beamerscope/MATLAB/'; % MAC
%root_path = 'C:\Users\Bene\Dropbox\Dokumente\OpticalDesign\Beamerscope\MATLAB\';

% rmpath(strcat(root_path, 'Matlab\TCC'))
addpath(genpath(strcat(root_path, 'ILLOPT_MATLABANDROID')))
addpath(genpath(root_path))
cd(strcat(root_path, 'ILLOPT_MATLABANDROID'))




%% System parameters
%% Set parameters of the acquisition process (in �m)
wavelength=0.530;       % systems center-wavelength
NAo = 0.25;             % aperture of microscopes objective
NAc = 0.5;              % maximum condenser aperture

pixel_eff = 0.6;%1.2;        % representing effective pixelsize of the camera
pixel_n = 61;           % number of simulated pixel

n_immers = 1.10;        % refractive index immmersion medium
n_obj = 1.20;           % refractive index object

%% Set parameters of the optimization process
sourceshape = 'MultiSegment';    % base for source shape optimization ('Zernike', 'MultiSegment', 'Archels')
n_coefficients = 48;             % number of zernikes used for optimizaton



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Programm starts here %%%%%%%%%%%%%%
%% Calculate spatial/frequency grid

simRadius=((pixel_n-1)*pixel_eff)/2;      % radius of simulated region in microns
xsim=-simRadius:pixel_eff:simRadius;      % create grid for simulation
zsim= 0;                                  % simulated z-position (not yet implemented)
usim=0;                                   %

kernel_limit = 50;


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

if(true)
    %% precalculate test-TCC and eigenvalues for desired parameterset
    [ IOsys ] = precalculatetcc( IOparams, xsim, zsim );
    
    %% assign eigenvalues/eigenfunctions
    titlestr = strcat('./result/', datestr(now, 'dd-mm,HH-MM-SS'),'_', IOparams.sourceshape, '_NAc', num2str(IOparams.NAc), '_NAo', '_rotAng', '_nCoefficients', num2str(IOparams.n_coefficients), '_xSize', num2str(size(xsim,2)));
    eigenvalue = IOsys.eigenvalue;
    %save(strcat(titlestr, 'eigenvalue.mat'), 'eigenvalue');
    eigenfunction = IOsys.eigenfunction;
    
    
    %% Optimization parameters
    costfun = 'MSE'; % NAE, MSE, STD, AD, CCP, AMM, IMMSE, NRMSE, SC, NCC, MSE, CPP
    opt_method = 'swarm';        % optimization algorithm to optimize soource shape ('swarm', 'grid')
    kernel_limit = 50;          % how many svd-kernels have to be tacken into account for simulation
    ill_method = sourceshape;         % which illumination source of choice?
    n_coefficients = 48;
    MaxTime = 100;
    kernelOrder = 1;
    ub = 1;
    lb = 0;
    StallIterLimit = 20;
    MaxIter = 40;
    kernel_limit = 50;
    debug = true;
    
    %% Setup initial coefficients
    x0 = rand(1, n_coefficients);
    x0 = ones(size(x0));
    
    %% Manual Mask
    % just in case you want to enable/disable specific patterns (i.e. only
    % outer circle when using
    %xmask= [1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 7];
    xmask = [1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    xmask = ones(size(x0));
    
    
    
    %% setup optimization parameters
    disp('------Optimize Illumination!---------------')
    OPTparams.x0 = x0; % initial optimization parameter set
    OPTparams.xmask = xmask;
    OPTparams.optmethod =  opt_method;  % searchmethod
    OPTparams.kernelOrder = kernelOrder;
    OPTparams.MaxTime = MaxTime;
    OPTparams.n_coefficients = n_coefficients;
    OPTparams.ub = ub;
    OPTparams.lb = lb;
    OPTparams.StallIterLimit = StallIterLimit;
    OPTparams.MaxIter = MaxIter;          % max. iteration to find optimum
    OPTparams.StallIterLimit = StallIterLimit;
    OPTparams.qualitymethod = costfun;
    
end


% load the images from the preprocessed database
cd('./resultOptNN/')
load('PreProcessedDataNN_filtered');
nFiles = size(complxObject, 3);

clear obj_complx_1d
clear obj_complx_2d
index = 1;

newsize = 128;
obj_complx_2d = zeros(nFiles*12, newsize, newsize, 2);
obj_complx_1d = zeros(nFiles*12, newsize*newsize*2);
x_opt1D = zeros(nFiles*12, 48);

for iImages=2%(1:nFiles)
    
    disp( [num2str(iImages), '/', num2str(nFiles)] )
    object = complxObject_filtered(:,:,iImages);
    
    %% important for intensities of segments; values between 0..1
    %% optimize illumination source
    [ x_opt,fval,output ] = optillsourceAndroid(ill_method, eigenvalue, eigenfunction, OPTparams, object);
    
    object = imresize(object, [newsize, newsize]);
    
    
    
    %% save data for later use
%     x_opt1D(:,iImages) = x_opt-0.5; %zerocenter the data for the NN
%     x_opt2D{iImages} = reshape(x_opt, [], 4)';
%     fval_i{iImages} = fval;
%     output_i{iImages} = output;
    
    
    
    
    % rotate the images and optimized illpatterns respectivly; this will
    % inrease the dataset size without recalculating the entire dataset
    % again; reason for that is the FT-theorem.
    for iAngles = 0:11
        
        
        %x_opt
        % rotate optimized pattern as well and bring back into 1D-vector
        % shape
        x_opt1D(index,:) = reshape(circshift(reshape(x_opt, [], 4)', [0 iAngles-1 ])', 1, []);
        
        % rotate object => generates more data!
        object_rotate = imrotate(object, iAngles*22.5, 'crop');
        
        if(true)
            % in this case use the real/imaginary part of the image as two
            % layer-image
            
            %obj_complx_2D(index,:,:,:) = cat(3, abs(object_rotate), angle(object_rotate));
            
            obj_complx_2d(index, :, :, :) = cat(3, abs(object_rotate), angle(object_rotate));
            
            
            % transform 2d map into 1d vector for NN
            % obj_magn_1d = reshape(obj_magn_2d, 1, []);
            % obj_phase_1d = reshape(obj_phase_2d, 1, []);
            
            %obj_abs_1d = reshape(obj_abs_2d, 1, []);
            %obj_angle_1d = reshape(obj_angle_2d, 1, []);
            %obj_complx_1d(index, :) = [obj_abs_1d obj_angle_1d];
            
            
            %obj_complx_1d(index, :) = [obj_abs_1d obj_angle_1d];
            
        else
            % in this case use the real/imaginary part of the image 
            % spectrum as two layer-image
            
            
            obj_ft = fftshift(fft2(object_rotate));
            obj_ft = imresize(obj_ft, [64, 64]);
            
            % decrease resolution of ft-spectrum
            % obj_ft = imresize(obj_ft, 0.25);
            
            % obj_magn_2d = abs(obj_ft);
            % obj_phase_2d = angle(obj_ft);
            
            obj_real_2d = real(obj_ft);
            obj_imag_2d = imag(obj_ft);
            
            
            % transform 2d map into 1d vector for NN
            % obj_magn_1d = reshape(obj_magn_2d, 1, []);
            % obj_phase_1d = reshape(obj_phase_2d, 1, []);
            
            obj_real_1d = reshape(obj_real_2d, 1, []);
            obj_imag_1d = reshape(obj_imag_2d, 1, []);
            
            
            % create complex variable for 1d object
            % obj_complx_1d(index, :) = [obj_magn_1d obj_phase_1d];
            obj_complx_1d(index, :) = [obj_real_1d obj_imag_1d];
        end
        
        
        if(false)
            %% generate illumination source for debugging with optimized pattern
            IOparams.Zcoefficient = x_opt1D(index,:);
            AERIALsys=tcc(xsim,zsim);
            AERIALsys.computesys(IOparams.sourceshape, IOparams);
            % imwrite(AERIALsys.Ic, strcat(currentfilename, num2str(iAngles), '.jpg'));
            
            % Calculate TCC by stacking the pupil overlaps in a 2D NxM^2 Matrix
            [TCC2D] = calculateTCC(AERIALsys);
            
            % Decomposing the TCC into its eigenvectors and eigenvalues
            [ef ev] = calculateSVD(AERIALsys);
            
            % Calculate the aerial image from the TCC and object-spectrum
            [IAerial] = (calculateImage(AERIALsys, object));
            
            
            %% generate illumination source for debugging with initial pattern
            IOparams.Zcoefficient = x0;
            AERIALsys=tcc(xsim,zsim);
            AERIALsys.computesys(IOparams.sourceshape, IOparams);
            % imwrite(AERIALsys.Ic, strcat(currentfilename, num2str(iAngles), '_init.jpg'));
            
            % Calculate TCC by stacking the pupil overlaps in a 2D NxM^2 Matrix
            [TCC2D] = calculateTCC(AERIALsys);
            
            % Decomposing the TCC into its eigenvectors and eigenvalues
            [ef ev] = calculateSVD(AERIALsys);
            
            % Calculate the aerial image from the TCC and object-spectrum
            [IAerial_init] = (calculateImage(AERIALsys, object));
            
            figure
            subplot(121)
            imagesc(IAerial_init)
            colormap gray
            axis square
            colorbar
            title('Init')
            
            
            subplot(122)
            imagesc(IAerial)
            colormap gray
            axis square
            colorbar
            title('Opt')
        end
        
        index = index + 1;
        
    end
    
    
end

% save('result_Dataset','-v7.3')


% for num_images = 1:nfiles
%     for num_rotation = 1:12
%
%     end
% end

% make sure data is normalized and zero-centered!



% prepare for TENSORFLOW - sorry for this inconvenience
nninputs = ((obj_complx_2d));
nnoutputs = x_opt1D;

% save data before processing it again
save('nninputs_realimag', 'nninputs', '-v7.3')
save('nnoutputs', 'nnoutputs', '-v7.3')

return


% try
%     nnoutputs = x_opt1D-.5; % zero-centre data for NN application (in case of Segments it's simply -0.5)
%     nnoutputs = ((nnoutputs));
% catch
% end


%% take care of the "wrong" data
% for spectrum
% nninputs_filtered = nninputs(~(sum(abs(nnoutputs), 2)>=23),:);
% nnoutputs_filtered = nnoutputs(~(sum(abs(nnoutputs), 2)>=23),:);

% for 2d real space


% for real space
in_max = max(nninputs, [], 2);
nninputs_filtered = nninputs(in_max>1e-1, :);
nnoutputs_filtered = nnoutputs(in_max>1e-1,:);

save('nninputs_realimag', 'nninputs_filtered', '-v7.3')
save('nninputs_realimag_2D', 'obj_complx_2d', '-v7.3')
save('nnoutputs', 'nnoutputs_filtered', '-v7.3')

% %% generate test data
%
% % get 10 percent from dataset as testdata
% percent_training = 10;
% size_dataset = size(nninputs,2);
%
% % generate Testdata
% index_i = randi([1, size_dataset],1,round(.1*size_dataset));
%
% nnoutput_test = nninputs(:, index_i);
% nnoutput_test =  nnoutputs(:, index_i);
%
% % save training data
% save('nninput_test', 'nninput_test')
% save('nnoutput_test', 'nnoutput_test')


return

%% do the inverse check:
%magnitude-part
imagesc(reshape(nninputs(1, 1:128^2), [128, 128]))
%phase-part
imagesc(reshape(nninputs(1, 128^2+1:end), [128, 128]))

% display original phase object
angle(ift(reshape(nninputs(1, 1:128^2), [128, 128]).*exp(1i*reshape(nninputs(1, 128^2+1:end), [128, 128]))))


%% check for weird values
imagesc(sqrt(abs(nninputs)))

%% save only phase
imagesc(nninputs(:,128^2/2+1:end))

%% save only magnitude
imagesc(nninputs(:,1:128^2/2))


