function [weights1D weightmap_ft] = getCircSpec(object_ft, n_segments, n_rings)

% this function takes an input spectrum and generates a circular
% segmentation map of the spectrum
%
% the vector holds its weights in n_rings and n_segments in each ring;
% the 2D matrix holds the map for debug purposes

segment_stack = generateSegments( size(object_ft), n_segments, n_rings);
weights1D = simplifySpectrum(object_ft, segment_stack);
      
weightmap_ft = zeros(size(object_ft));  
% plot simplified object spectrum
for i = 1:(n_segments*n_rings)
    weightmap_ft = weightmap_ft + weights1D(i).*segment_stack(:,:,i);
end

% %% debug 
% figure
% subplot(1,2,1)
% imagesc(log(object_ft))
% colormap gray, axis square
% 
% subplot(1,2,2)
% imagesc(log(weightmap_ft))
% colormap gray, axis square