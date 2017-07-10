function TCC = integrateSOCSPAR(Po, J, S)
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


disp('integrate Pupil');

xSize = size(J, 2).^2;
ySize = size(J, 1).^2;

% some stuff for shifting the pupils
N = size(J, 2);
midway=(N+1)/2;

% generate Illumination aperture
TCC = zeros(xSize, ySize);

% get integration boundaries
Po_binary = im2bw(abs(Po)+1,1);  % Make sure, that even if the pupil is complex or negative, the entire radius is determined
r_NAo = round(sum(Po_binary(:,round(size(Po_binary,1)/2)))/2)+1;    % radius of objectives aperture in pixel
r_NAc = r_NAo * S;                                                  % radius of condensers aperture in pixel


debug = true;




%%%%%%Calculate transmission cross coefficients%%%%%%
parfor x=1:xSize
    for y=1:ySize
        
        P_1=Po;
        P_2=conj(Po);
        
        index_m=mod(x-1,N)+1-midway;
        index_n=floor((x-1)/N)+1-midway;
        index_p=mod(y-1,N)+1-midway;
        index_q=floor((y-1)/N)+1-midway;
        
        
        %check if the both objecitve puils have intersections
%         if ((sqrt(index_m^2+index_n^2)>(r_NAo+r_NAc)/2) | (sqrt(index_p^2+index_q^2)>(r_NAo+r_NAc)/2))
%             TCC(x,y) = 0;
%        %elseif(sqrt((shiftmatrix(center_m)-shiftmatrix(center_p))^2+(shiftmatrix(center_n)-shiftmatrix(center_q))^2)>(2*r_NAo))             TCC(x,y) = 0;
            if(false)
                
            else
            
            if (index_m>0)
                P_1=[Po(index_m+1:N,1:N);zeros(index_m,N)];
            else
                P_1=[zeros(abs(index_m),N);Po(1:N+index_m,1:N)];
            end
            
            if (index_n>0)
                P_1=[P_1(1:N,index_n+1:N),zeros(N,index_n)];
            else
                P_1=[zeros(N,abs(index_n)),P_1(1:N,1:N+index_n)];
            end
            
            
            if (index_p>0)
                P_2=[Po(index_p+1:N,1:N);zeros(index_p,N)];
            else
                P_2=[zeros(abs(index_p),N);Po(1:N+index_p,1:N)];
            end
            if (index_q>0)
                P_2=[P_2(1:N,index_q+1:N),zeros(N,index_q)];
            else
                P_2=[zeros(N,abs(index_q)),P_2(1:N,1:N+index_q)];
            end
            TCC(x,y)=sum(sum(J.*P_1.*P_2));
        end
        
    end
end



