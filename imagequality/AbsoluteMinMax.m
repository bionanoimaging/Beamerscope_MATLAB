function [ AMM ] = AbsoluteMinMax( distImg )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
I_min = min(min(distImg));
I_max = max(max(distImg));
AMM = (I_max-I_min)./(I_max+I_min);

end

