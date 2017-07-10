function displayDPC( recparams )
%This script shows a specific aperture on the iphone screen; It uses the
%Psychtoolbox 3.0 to display the pattern full-screen on the second screen.
%The parameters set here are adequate for the current configuration. It
%uses a second monitor which is mirrored to the iphone (which itself is
%configured as a second display using the APP DesktopX from slpashtop) 
%
% Parameters like xpos/ypos have to be calibrated in advance!
% 
% Input: NAc - it represents the numerical aperture used by the condensor.
% Standard is NAc = 0.9.
%% 

% LCD pixel number
dim_desktop = [1920 1200];  % resolution of the 2nd monitor
dim_iphone = [960 640];     % resolution of the iphone-LCD at 326 dpi

% Central position of circular aperture; callibration by centering the
% aperuter while observing the BFP through the bertrand lens
xpos = 980;
ypos = 610;

% Radius of Condenser NA in Pixel == 0.9NA (in 2nd Screen coordinates!)
radius = 392.5*(recparams.NAc/0.9);

% Angles of semi circles
angle_semi = [0 pi pi*3/2 1/2*pi];    % order: left right, top bottom

%% Script starts here
x = linspace(1, dim_desktop(1), dim_desktop(1));
y = linspace(1, dim_desktop(2), dim_desktop(2));

[xx yy] = meshgrid(x',y');

% get shape of aperture
%theImage = circular_pupil( xx, yy, xpos, ypos, radius );

for(i=1:4)
    theImage(:,:,i) = semicircular_pupil( xx, yy, xpos, ypos, radius, angle_semi(i) );
end


% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');
Screen('Preference','SkipSyncTests', 1);
Screen('Preference', 'VisualDebugLevel', 0);

ShowCursor;
% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');



for(j=0:size(theImage,3)-1)
    
    i = mod(j, 4)+1;
    circle = imadjust(double(theImage(:,:,i)));
    
    % Get the size of the image
    [s1, s2, s3] = size(circle);
    
    % Here we check if the image is too big to fit on the screen and abort if
    % it is. See ImageRescaleDemo to see how to rescale an image.
    if s1 > screenYpixels || s2 > screenXpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        return;
    end

    % Make the image into a texture
    imageTexture = Screen('MakeTexture', window, circle);
    
    % Draw the image to the screen, unless otherwise specified PTB will draw
    % the texture full size in the center of the screen. We first draw the
    % image in its correct orientation.
    Screen('DrawTexture', window, imageTexture, [], [], 0);
    
    % Flip to the screen
    Screen('Flip', window);
    
    % Wait for two seconds
    WaitSecs(2);
end

sca;