function dispPattern( dispTime, dispMethod, dispParams)
%UNTITLED3 Summary of this function goes here
%   dispTime - time to display Pattern
%   dispMethod - pattern type (Zernike, Multisgment, Archels)
%   dispParams - Coefficients for dispMethod


% Display Illumination shape

% LCD pixel number
dim_desktop = [1920 1200];
dim_iphone = [960 640];     %326 dpi
% Radius of Condenser NA in Pixel == 0.9NA
radius_NA = 192.5;
NA_physical = 0.9;
NA_used = 0.5;
radius = (radius_NA/NA_physical)*NA_used;

NAo=0.5;
NAc=NA_used;
% load display calibration-parameter
Data()

% create grid
x = linspace(1, dim_desktop(1), dim_desktop(1));
y = linspace(1, dim_desktop(2), dim_desktop(2));
[xx yy] = meshgrid(x',y');
circle = circular_pupil( xx, yy, xpos, ypos, radius);

if(strcmp(dispMethod, 'Multisegment'))
    pattern = (SegCalc( circle, dispParams, NAo, NAc ));
elseif(strcmp(dispMethod, 'Zernike'))
    pattern = (SegCalc( circle, dispParams ));
elseif(strcmp(dispMethod, 'Archels'))
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
% close all

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


% Make the image into a texture
imageTexture = Screen('MakeTexture', window, pattern);

% Draw the image to the screen, unless otherwise specified PTB will draw
% the texture full size in the center of the screen. We first draw the
% image in its correct orientation.
Screen('DrawTexture', window, imageTexture, [], [], 0);

% Flip to the screen
Screen('Flip', window);

% Wait for two seconds
WaitSecs(dispTime);


% Clear the screen
sca;
end

