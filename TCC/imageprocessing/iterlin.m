function it = iterlin(rng)
it_type = size(rng, 2);
if it_type == 2
    rng = rng(:, 2) - rng(:, 1) + 1;
elseif it_type > 2
    error('bad iterand format.');    
end
it = 1 : prod(rng);