function [ TCC4D ] = reshapeTCC( TCC2D )
%reshapeTCC it simply reshapes the 2D TCC into the appropriate 4D shape
%   a 2D optical system follows in a 4D Transferfunction; For faster
%   computation the TCC comes as a 2D matrix; For better assignment it has
%   to be "reshaped" 
%   Taken from here http://stackoverflow.com/questions/20336288/for-loop-to-split-matrix-to-equal-sized-sub-matrices

%   convert 2D-TCC -> 4D-TCC

N_mask = sqrt(size(TCC2D, 1)); 
TCC4D = permute(reshape(TCC2D, size(TCC2D, 1), N_mask, []), [2 1 3]);
TCC4D = permute(reshape(TCC4D, N_mask, N_mask, [], size(TCC4D, 3)), [2 1 3 4]);

end

