function size = vnorm(vec,dim)

% Basic 2-norm of a vector or row/columns of a matrix.

size = sqrt(sum(vec.^2,dim));