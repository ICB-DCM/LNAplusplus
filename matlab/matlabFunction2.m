function f=matlabFunction2(m, varargin)
% this function takes a symbolic matrix and converts it to a matrix of
% symbols, then calls matlabFunction.  this is for some reason much faster
% than directly calling matlabFunction on the symbolic matrix
% returns a function handle

tmp = reshape( m(:), size(m));
f = matlabFunction(tmp, varargin{:});

end