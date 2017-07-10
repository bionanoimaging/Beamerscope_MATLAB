function displayPAT( recmethod, recparams )
%displayPAT Simply switches between different illumination configurations
% %   
% recmethod = 'DPC', 'FPM'
% recparams = struct(NA_c,..)

if  strcmp(recmethod,'DPC')
    displayDPC(recparams)
elseif strcmp(recmethod,'FPM')
    displayFPM(recparams)
end

