function [ TCC2D ] = reshapeTCC4d2d( TCC4D )
%reshapeTCC it simply reshapes the 2D TCC into the appropriate 4D shape
%   a 2D optical system follows in a 4D Transferfunction; For faster
%   computation the TCC comes as a 2D matrix; For better assignment it has
%   to be "reshaped" 
%   Taken from here http://stackoverflow.com/questions/20336288/for-loop-to-split-matrix-to-equal-sized-sub-matrices

%   convert 4D-TCC -> 2D-TCC


[s1,ign1,s3,ign2] = size(TCC4D);
TCC2D = reshape(permute(TCC4D,[1 3 2 4]),s1*s3,[]);

end

