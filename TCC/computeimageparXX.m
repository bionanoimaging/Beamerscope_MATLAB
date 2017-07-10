function aerial_fre=computeimageparXX(TCC, obj)
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


% Inspired by the work of "Computational Lithography" XU MA AND GONZALO R. ARCE
%   "PSM_svd.m" Two-phase PSM optimization using the singular value decomposition
%   model in partially coherent imaging system.

% for further info see here: http://de.mathworks.com/matlabcentral/answers/261983-sliced-variables-in-parfor-loop-any-chance-to-optimize

disp(strcat('Start computing the image using matrix optimization. Please wait...'));

%% mathwork solution
tic

N_mask = size(obj,2);
pz_fre=(fftshift(fft2(obj)));

[index_x, index_y] = ndgrid(1:N_mask^2);

[index_1,index_2] =  ind2sub([N_mask,N_mask],  index_x);
[index_3,index_4] =  ind2sub([N_mask,N_mask],  index_y);

index_m = (mod(index_1-index_3,N_mask) + 1).';
index_n = (mod(index_2-index_4,N_mask) + 1).';

vals=bsxfun(@times, TCC,pz_fre(:));
vals=bsxfun(@times,vals,pz_fre(:)');

idx=(vals~=0);

subs=[index_m(idx), index_n(idx)];
vals=vals(idx);

aerial_fre = accumarray(subs,vals,[N_mask,N_mask]);
    
    
toc



end
