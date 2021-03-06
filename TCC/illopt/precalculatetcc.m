function [ IOsys ] = precalculatetcc( IOparams, xsim, zsim )
% generates the TCC for a given microscope configuration 
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.
%
% function [ IOsys ] = precalculatetcc( IOparams, xsim, zsim )
% inputs: IOparams - Parameters for the microscope setup 
%       xsim - simulation grid in x/y-direction
%       zsim - simulation grid in z-direction
% outputs:
%       IOsys - generated microscope system with physical pupil shapes and
%       its corresponding eigenfunctions/eigenvalues




% if ZERNIKE
if(strcmp(IOparams.sourceshape, 'Zernike'))
    
    Zernike_coefficients = eye(IOparams.n_coefficients);
    
    % allocate memory
    eigenfunction = zeros(size(IOparams.xsim,2), size(IOparams.xsim,2), IOparams.kernel_limit, IOparams.n_coefficients);
    eigenvalue = zeros(IOparams.kernel_limit, IOparams.n_coefficients);
    
    for(numiter = 1:IOparams.n_coefficients)%size(Zernike_coefficients,1))
        
        % select the i-th Zernike coefficient
        IOparams.Zcoefficient = Zernike_coefficients(numiter,:);
        
        % allocate memory for i-th coefficient and calculate system
        % (apertures)
        IOsys=tcc(xsim,zsim);
        IOsys.computesys(IOparams.sourceshape, IOparams);
        
        %% Major Calculation is done here
        % Calculate TCC by stacking the pupil overlaps in a 2D NxM^2 Matrix
        [TCC2D] = calculateTCC(IOsys);
        
        % Decomposing the TCC into its eigenvectors and eigenvalues
        IOsys.N = 50; % limit the number of kernels generated by the svd to a definite number, otherwise the process of archiving the database wont work
        [ef ev] = calculateSVD(IOsys);
        
        % save all eigenfunctions in an array for later use
        eigenfunction(:,:,:,numiter) = ef;
        eigenvalue(:,numiter) = ev;
        
        
    end
    
    
elseif(strcmp(IOparams.sourceshape, 'MultiSegment'))
    
    Segment_Intensities = eye(IOparams.n_coefficients);
    
    % allocate memory
    eigenfunction = zeros(IOparams.pixel_n, IOparams.pixel_n, IOparams.kernel_limit, IOparams.n_coefficients);
    eigenvalue = zeros(IOparams.kernel_limit, IOparams.n_coefficients);
    
    disp('Start to compute the TCC in advance...');
    for(numiter = 1:IOparams.n_coefficients)%size(Zernike_coefficients,1))
        disp(strcat('TCC: ', num2str(numiter), '/', num2str(IOparams.n_coefficients)));
        IOparams.Zcoefficient = Segment_Intensities(numiter,:);
        
        IOsys=tcc(xsim,zsim);
        IOsys.computesys(IOparams.sourceshape, IOparams);
        
        %% Major Calculation is done here
        % Calculate TCC by stacking the pupil overlaps in a 2D NxM^2 Matrix
        [TCC2D] = calculateTCC(IOsys);
        
        IOsys.N = 50;
        % Decomposing the TCC into its eigenvectors and eigenvalues
        [ef ev] = calculateSVD(IOsys);
        
        eigenfunction(:,:,:,numiter) = ef;
        eigenvalue(:,numiter) = ev;
        
    end
    
    

   
end

IOsys.eigenfunction = eigenfunction;
IOsys.eigenvalue = eigenvalue;
end

