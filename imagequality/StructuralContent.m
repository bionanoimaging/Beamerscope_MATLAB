%Program for Structural Content Calculation

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

function SC = StructuralContent(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

SC = sum(sum(origImg .* origImg)) / sum(sum(distImg .* distImg));