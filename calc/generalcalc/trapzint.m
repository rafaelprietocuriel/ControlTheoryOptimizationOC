function z = trapzint(x,y,dim)
%TRAPZINT  Trapezoidal numerical integration.
%
%   Z = TRAPZINT(X,Y,DIM) computes the integral of Y with respect to the vector of intervals X using
%   the trapezoidal method.  Length of vector X must be one element less then the length of DIM dimention of array Y
%   length(X) = size(X,DIM) - 1 .  TRAPZINT operates along this
%   dimension.
%
%   Z = TRAPZINT(X,Y) integrates across the first dimension DIM = 1
%   of Y. The length of X must be the same as size(Y,DIM) - 1.
%
%   For N-D arrays, TRAPZINT(Y) works across the first dimension.
%
%   Make sure x and y are column vectors, or y is a matrix.
%
% Anton O. Belyakov 
% Revision 2010/03/05

switch nargin
  case 3
    nd = ndims(y);
    if dim > nd, error('Too high index Dim')  
    elseif dim > 1 % put the dimention dim in the first place 
      y = permute(y,[dim 1:dim-1 dim+1:nd]);
    elseif dim < 1, error('Index must be Dim > 0')
    end   
  case 2
    if isscalar(y)
      if isempty(x)
        z = 0;
        return
      end
      z = trapzsum(x,y);
      return
    else
      nd = ndims(y);
    end
  case 1
    z = trapzsum(x);
    return
  otherwise
    if nargin > 3 
      error('Too many input arguments')
    else 
      error('Not enough input arguments')
    end
end
 
m = size(y,1);
if length(x) ~= m-1, error('length(x) must equal size(y,dim) - 1.'); end
if m == 0, error('y must be nonempty'); end

if nd == 2 % number of dimentions (not less than 2, see help ndim) 
  if size(y,2) == 1, z = trapzint1q(x,y);
  else z = trapzint2q(x,y);
  end
else z = trapzint3q(x,y);   
end
