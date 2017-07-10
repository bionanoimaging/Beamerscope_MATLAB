%% Generate a dataset with complex objects and its corresponding optimized
% illumination shapes using the TCC
% 
%
% Scipt which reads all images and pre-processes them to have same
% properties and damped edges
%
% The software is for private use only and gives no guarrantee, that it's
% working as it should! 
% 
%
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.



%% read all images inside the dataset
imfolder = 'C:\Users\Bene\Downloads\DatasetMLcrop\'; % WINDOWS
imfolder = '\\Q-SERVER\USBDisk1\DatasetMLcrop\'; % WINDOWS SERVER
imfolder = '/Users/Bene/Dropbox/Dokumente/OpticalDesign/Beamerscope/DatasetML/'; % MAC


% WINDOWS
root_path = 'C:\Users\Bene\Dropbox\Dokumente\OpticalDesign\Beamerscope\MATLAB\';

% MAC
root_path = '/Users/Bene/Dropbox/Dokumente/OpticalDesign/Beamerscope/MATLAB/';

cd(strcat(root_path, 'ILLOPT_MATLABANDROID/resultOptNN'));

% look for all images in folder
imagefilesJPG = dir(strcat(imfolder, '*.jpg'));
imagefilesPNG = dir(strcat(imfolder, '*.png'));
nFiles = length(imagefilesJPG);    % Number of files found

% parameters for acquisition 
lambda = 550E-9;    % system centre-wavelength in m
n_immers = 1.2;     % refractive index immersion media
n_object = 1.22;     % refractive index object
objThickness = 10E-6 % thickness of the object

finalSizeImage = 128;    % All images must have this size
scaling_factor = 4; % the spectrom should be scaled DOWN by a factor of 4 => mxn=32x32, reduces noise and computational time

%% generate Input and Output Vector for Neural network
index = 1;

% initliaze variable
complxObject = zeros(finalSizeImage, finalSizeImage, nFiles);

%% preprocess all images
for i=(1:nFiles)
    disp(strcat(num2str(i), '/', num2str(nFiles)));
    
    % read current image
    currentfilename = strcat(imfolder, imagefilesJPG(i).name);
    
    try
        currentImage = imread(currentfilename);
        
        % convert image to double
        try
            currentImage = im2double(rgb2gray(currentImage));
        catch
            currentImage = im2double(currentImage);
        end
        
        % normalize image to 0..1
        currentImage = currentImage - min(min(currentImage));
        currentImage = currentImage./max(max(currentImage));
        
        % make sure that image is square
        sizeImage = min(size(currentImage));
        
        % handle images which are larger or smaller than the
        if finalSizeImage < sizeImage % if bigger, extract centre
            if finalSizeImage < 0.75*sizeImage % if 2*bigger, first shrink
                currentImage = imresize(currentImage, finalSizeImage/sizeImage*1.25);
            end
            currentImage = extract(currentImage, [finalSizeImage, finalSizeImage]);
        elseif finalSizeImage > sizeImage % if smaller
            currentImage = imresize(currentImage, finalSizeImage/sizeImage);
            currentImage = extract(currentImage, [finalSizeImage, finalSizeImage]);
        end
        
        % damp edges before - better results when doing the FT()
        objectProfile = currentImage;% DampEdge(currentImage,0.05,2,1,2);
        % cat(3, ft(objectProfile), ft(currentImage));
        
        
        % convert gray-values into phase-information
        oplObject = (2*pi/lambda)*(n_object-n_immers)*objectProfile*objThickness;
        
        % convert gray-values into phase-information
        complxObject(:,:,i) = 0.01*objectProfile*exp(1i*oplObject);
        
    catch
        print 'wrong image format'
    end
    
    
    
end

% filter due to wrong pixelvalues 
complxObject_filtered = complxObject(:,:,max(max(angle(complxObject), [], 2), [], 1)>0);

save('PreProcessedDataNN', 'complxObject');
save('PreProcessedDataNN_filtered', 'complxObject_filtered');
