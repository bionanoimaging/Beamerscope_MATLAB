function aerial_fre=computeimage(self, TCC, obj, params, computeDevice)
% COMUTETCC: Computes images from the imaging system's transmission cross
% coefficient following Hopkins method
% functions contained in the microlith object and the specimen properties.
%   tcc=computetcc(MLOBJ, computeDevice) computes the TCC
%   The pupil of the objective and its conjugate will be shifted accross
%   the plane of normalized optical coordinates. The intersection of the
%   two and the effictive illumination source of this setup is summarized
%   and assigned to the X/Y-pixel. The Result is a 4D-matrix and can be
%   used to simulate image formation i
%
%   computeDevice = 'CPU', 'multiCPU', or 'GPU'.
%   'multiCPU' option computes images across the focus in parallel if the
%   MATLAB parallel computing toolbox is installed and a matlabpool is
%   available.
%
%
%   Written by Shalin Mehta, www.mshalin.com
%   License: GPL v3 or later.
%   Modified by Benedict Diederich

% This file is part of microlith package.
%
%     microlith is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     microlith is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with microlith.  If not, see <http://www.gnu.org/licenses/>.

if(~isa(self,'tcc'))
    error('This function requires an object of the microlith class.');
end


% Inspired by the work of "Computational Lithography" XU MA AND GONZALO R. ARCE
%   "PSM_svd.m" Two-phase PSM optimization using the singular value decomposition
%   model in partially coherent imaging system.

disp(strcat('Start computing the image. Please wait...'));


startTime = cputime;

tic

N_obj = size(obj,2);

index_x = repmat(1:N_obj^2,[N_obj^2 1]);
index_x = index_x(:);
index_y = repmat(1:N_obj^2, 1, N_obj^2)';

index_1=mod(index_x-1,N_obj)+1;
index_2=floor((index_x-1)/N_obj)+1;
index_3=mod(index_y-1,N_obj)+1;
index_4=floor((index_y-1)/N_obj)+1;

index_m = (mod(index_1-index_3,N_obj) + 1)';
index_n = (mod(index_2-index_4,N_obj) + 1)';

aerial=zeros(N_obj,N_obj);
aerial_fre=zeros(N_obj,N_obj);
obj_fre=(fftshift(fft2(obj)));

for idx=1:N_obj^4    
    aerial_fre(index_m(idx),index_n(idx))=aerial_fre(index_m(idx),index_n(idx))+...
        TCC(index_x(idx),index_y(idx))*(obj_fre(index_1(idx),index_2(idx)))*...
        conj(obj_fre(index_3(idx),index_4(idx)));
end  

toc




end
