function TCC = integratePupilPar(Po, J, S)
% integratePupilPar implements the discretesized integration of the
% transmission cross coeficient. It calculates the sum of the overlapping
% areas from the source and objective pupil depending on the object's
% spectrum frequencies
%
%   Written by Shalin Mehta, www.mshalin.com
%   License: GPL v3 or later.
%   Modified by Benedict Diederich

% I assume that all the matrices are defined over the same frequency grid.
% This function is not immune to errors. Error checks must be
% implemented by the parent function. The error checks are ommited
% so that the function can be run on GPU without having to transfer any
% data to CPU.


% The spectral support (or spatial cutoff) of simulation must be able
% to accommodate the largest spatial frequency that is captured.
% Therefore, there is no need to pad either objspecG or Po to find the
% product of shifted objspecG and Po.

%Assignment in this manner ensures that new variables are on GPU if Ns
%is on GPU. Note that Ns is real, whereas objSpectrum and Po are
%complex.


%waith=waitbar(loopvar/Ns,['Computing image intensity by summing over the source:' num2str(100*loopvar/Ns,2) '%']);


% For each source point, calculate what indices from the unshifted
% spectrum map to what indices of shifted spectrum.
%---------------------------------------

% Where do shifted and unshifted spectra start given the sign of mshift.
%---------------------------------------------------------------

%     % Above is the vectorized form of the following logic for each mshift.
%     % It is vectorized to run efficiently on GPU.
%          if(mshift<0)
%              mSpectrumidx= (-mshift+1):L;
%              mShiftedSpectrumidx= 1:(L+mshift);
%
%          else
%              mSpectrumidx=1:(Len-mshift);
%              mShiftedSpectrumidx=(1+mshift):Len;
%          end

% Likewise for nshift
%--------------------------

% Need to vectorize the while loop. The idea is:
% Create a stacked shifted object spectrum with dimensions: (fx,fy,S).
% Multiply along (fx,fy) dimensions using Po with bsxfun.
% Obtain IFFT along (fx,fy) dimensions.
% Mag. square and sum along S.

disp('integrate Pupil');

xSize = size(J, 2);
ySize = size(J, 1);


% generate Illumination aperture
TCC = zeros(xSize, xSize, xSize, xSize);

% get integration boundaries
Po_binary = im2bw(abs(Po)+1,1);  % Make sure, that even if the pupil is complex or negative, the entire radius is determined
r_NAo = round(sum(Po_binary(:,round(size(Po_binary,1)/2)))/2)+1;    % radius of objectives aperture in pixel
r_NAc = r_NAo * S;                                                  % radius of condensers aperture in pixel

shiftmatrix = (-round(xSize/2)+1:round(xSize/2)-1);

debug = true;

parfor center_m = 1:xSize;         % shift pupil P_1 along X-frequency
    for center_p = 1:xSize;	% shift pupil P_2 along X-frequency
        for center_n = 1:ySize;             % shift pupil P_1 along Y-frequency
            for center_q = 1:ySize;     % shift pupil P_2 along Y-frequency
                
                P1=Po;
                P2=Po;
                
%                 %check if the both objecitve puils have intersections
%                 if ((sqrt(shiftmatrix(center_m)^2+shiftmatrix(center_n)^2)>(r_NAo+r_NAc)/2) | (sqrt(shiftmatrix(center_p)^2+shiftmatrix(center_q)^2)>(r_NAo+r_NAc)/2))
%                     TCC(center_m, center_p, center_n, center_q) = 0;
%                     elseif(sqrt((shiftmatrix(center_m)-shiftmatrix(center_p))^2+(shiftmatrix(center_n)-shiftmatrix(center_q))^2)>(2*r_NAo)))
%                        TCC(center_m, center_p, center_n, center_q) = 0;
                    if(false)
                else
                    
                    
                    
                    %decenter Pupil P2 in y-direction
                    if(shiftmatrix(center_q)>=0)
                        P2 = [Po(shiftmatrix(center_q)+1:xSize,1:xSize);zeros(shiftmatrix(center_q),xSize)];
                    end
                    if (shiftmatrix(center_q)<0)
                        P2 = [zeros(abs(shiftmatrix(center_q)),xSize);Po(1:xSize+shiftmatrix(center_q),1:xSize)];
                    end
                    
                    %decenter Pupil P2 in x-direction
                    if(shiftmatrix(center_p)>=0)
                        P2 = [P2(1:xSize,shiftmatrix(center_p)+1:xSize), zeros(xSize,shiftmatrix(center_p))];
                    end
                    if (shiftmatrix(center_p)<0)
                        P2=[zeros(xSize,abs(shiftmatrix(center_p))),P2(1:xSize,1:xSize+shiftmatrix(center_p))];
                    end
                    
                    
                    
                    
                    %decenter Pupil P1 in y-direction
                    if(shiftmatrix(center_n)>=0)
                        P1 = [Po(shiftmatrix(center_n)+1:xSize,1:xSize);zeros(shiftmatrix(center_n),xSize)];
                    end
                    if (shiftmatrix(center_n)<0)
                        P1 = [zeros(abs(shiftmatrix(center_n)),xSize);Po(1:xSize+shiftmatrix(center_n),1:xSize)];
                    end
                    
                    %decenter Pupil P1 in x-direction
                    if(shiftmatrix(center_m)>=0)
                        P1 = [P1(1:xSize,shiftmatrix(center_m)+1:xSize), zeros(xSize,shiftmatrix(center_m))];
                    end
                    if (shiftmatrix(center_m)<0)
                        P1=[zeros(xSize,abs(shiftmatrix(center_m))),P1(1:xSize,1:xSize+shiftmatrix(center_m))];
                    end
                    
                    P2 = conj(P2);
                    
                    
                    TCC(center_m, center_p, center_n, center_q) = sum(sum(P1.*P2.*J));
                    
                    
%                     if(debug)
%                         %%% Debugging
%                         
%                         figure(1)
%                         subplot(151)
%                         imagesc(P1)
%                         axis square
%                         colormap gray
%                         title 'P_1'
%                         
%                         subplot(152)
%                         imagesc(P2)
%                         axis square
%                         colormap gray
%                         title 'P_2'
%                         
%                         subplot(153)
%                         imagesc(J)
%                         axis square
%                         colormap gray
%                         title 'J'
%                         
%                         subplot(154)
%                         imagesc(J.*P1.*P2)
%                         axis square
%                         colormap gray
%                         title 'Overlap'
%                         
%                         if(mod(center_n+1,xSize)==0)
%                             subplot(155)
%                             imagesc(squeeze(TCC(center_m, :, center_p, :)))
%                             axis square
%                             colormap gray
%                             titletcc = strcat('TCC at m=', num2str(center_m), ' p=', num2str(center_p));
%                             title(titletcc)
%                         end
%                         
%                         drawnow
%                     end
                end                              
            end           
        end
        %disp(strcat(num2str(center_p), ' - ', num2str(center_q)));
    end
end





% Following is the normalization according to Martin's book. We do not use it, because the
% following leads to divide by zero for dark-field system.
%         normfactor=abs(Po).^2.*abs(IcG);
%         normfactor=sum(normfactor(:));
%         if(normfactor) %If normafactor == 0, we are imaging with dark-field system.
%             imgG=imgG/sum(normfactor(:));
%         end

end
