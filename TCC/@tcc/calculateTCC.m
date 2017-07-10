function [TCC M N] = calculateTCC(self)
% This function calculates the TCC according to the paper:
% Aerial Image Simulation for partial coherent system with programming 
% development in MATLAB
% Written by Benedict Diederich, benedict.diederich@leibniz-ipht.de
% www.nanoimaging.de
% License: GPL v3 or later.
%
%
% 
% input:
% - Source - 2D matrix with desired source pattern
% - p - 2D matrix containing the microscopes objective function
%
%
% output:
% - TCC - 2D TCC representation; each TCC is stacked into one row
% - M, N - Dimensions of the TCC

if(~isa(self,'tcc'))
    error('This function requires an object of the microlith class.');
end

M=size(self.Tfun,2);
N=numel(find(self.Ic~=0)); % number of point sources
[row,col]=find(self.Ic~=0); % position of point sources
TCC=zeros(N,M^2);
Ic_norm=self.Ic/sum(sum(self.Ic));


%disp('%%%%%%Calculate TCC%%%%%%%')
%tic
%parpool(2)
% To get the singular matrix P, stacking of pupil function has been done 
% which is shifted according to illumination

parfor i=1:N
CP=circshift(self.Tfun,round([row(i)-(M+1)/2,col(i)-(M+1)/2]));
TCC(i, :)=sqrt(Ic_norm(row(i),col(i)))*reshape(CP',1,[]); % singular matrix P
end
%toc

self.M = M;
self.N = N;
self.TCC = TCC;

end
