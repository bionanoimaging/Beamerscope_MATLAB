function Inorm = normminmax(I)
Inorm = I-min(min(I));
Inorm = Inorm./max(max(Inorm));
end
