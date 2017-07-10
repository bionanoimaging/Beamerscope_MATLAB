function TCC = integrateSOCS(P,J,S)
% integratePupilPar implements the discretesized integration of the
% transmission cross coeficient. It calculates the sum of the overlapping
% areas from the source and objective pupil depending on the object's
% spectrum frequencies
%
%   Written by Shalin Mehta, www.mshalin.com
%   License: GPL v3 or later.
%   Modified by Benedict Diederich

% P = Pupilfunction of the microscope objective [Mask x Mask]
% J = Pupilfunction of the lightsource [Mask x Mask]
% S = Coherence factor NAc/NAo

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

disp('integrate SOCS');

N = size(P, 2);
midway=round((N+1)/2);

% generate Illumination aperture
TCC = zeros(N^2, N^2);


startTime = cputime;

limit_x = N^2;
limit_y = N^2;

% get integration boundaries
Po_binary = im2bw(abs(P)+1,1);  % Make sure, that even if the pupil is complex or negative, the entire radius is determined
radius_Po = round(sum(Po_binary(:,round(size(Po_binary,1)/2-1)))/2);
radius_Pc = radius_Po * S;

if(matlabpool('size') < 2)
    matlabpool(4)
end

parfor x=1:limit_x
    P_1=P;
    P_2=conj(P);
    for y=1:limit_y
        
        % diffraction orders in X-direction
        index_m=mod(x-1,N)+1-midway;
        index_p=floor((x-1)/N)+1-midway;
        
        % diffraction orders in Y-direction
        index_n=mod(y-1,N)+1-midway;
        index_q=floor((y-1)/N)+1-midway;
        
        % First integration boundary condition: Po&Po^* have to have an overlap with Pc
        if( sqrt((index_m)^2+(index_p)^2) > (radius_Po + radius_Pc) || sqrt((index_n)^2+(index_q)^2) > (radius_Po + radius_Pc) )
            TCC(x,y)=0;
        elseif( sqrt((index_m-index_n)^2+(index_p-index_q)^2) > 2*radius_Po)
            TCC(x,y)=0;
        else
            
            
            if (index_m>0)
                P_1=[P(index_m+1:N,1:N);zeros(index_m,N)];
            else
                P_1=[zeros(abs(index_m),N);P(1:N+index_m,1:N)];
            end
            
            if (index_p>0)
                P_1=[P_1(1:N,index_p+1:N),zeros(N,index_p)];
            else
                P_1=[zeros(N,abs(index_p)),P_1(1:N,1:N+index_p)];
            end
            
            
            if (index_n>0)
                P_2=[P(index_n+1:N,1:N);zeros(index_n,N)];
            else
                P_2=[zeros(abs(index_n),N);P(1:N+index_n,1:N)];
            end
            if (index_q>0)
                P_2=[P_2(1:N,index_q+1:N),zeros(N,index_q)];
            else
                P_2=[zeros(N,abs(index_q)),P_2(1:N,1:N+index_q)];
            end
            TCC(x,y)=sum(sum(J.*P_1.*P_2));
%             
%             %%% Debugging
%             
%             figure(1)
%             subplot(141)
%             imagesc(P_1)
%             axis square
%             colormap gray
%             title 'P_1'
%             
%             subplot(142)
%             imagesc(P_2)
%             axis square
%             colormap gray
%             title 'P_2'
%             
%             subplot(143)
%             imagesc(J)
%             axis square
%             colormap gray
%             title 'J'
%             
%             subplot(144)
%             imagesc(J.*P_1.*P_2)
%             axis square
%             colormap gray
%             title 'Overlap'
%             
        end
        
        drawnow
    end
    
end
%disp(x);

%TCC = imrotate(TCC, 90);

endTime = cputime;

disp(strcat('Processing time is: ', num2str(startTime-endTime)));
end

