function varargout = iterget(rng, it)
assert(nargout == size(rng, 1), 'bad number of output iterands.');
it_type = size(rng, 2);
if it_type == 1
    dim = rng;
elseif it_type == 2
    dim = rng(:, 2) - rng(:, 1) + 1;
elseif it_type > 2
    error('bad iterand format');
end
[varargout{1 : nargout}] = ind2sub(dim', it);
if it_type == 2
    varargout = num2cell(cell2mat(varargout) + rng(:, 1)' - 1);
end