%toLinear
% This functino converts a symmetrix covariance matrix into a column-major-ordered
% vector for use in specifying the initial conditions for the LNA
% simulation
function f = toLinear( X )

f = X(find(triu(ones(size(X)))));

end

