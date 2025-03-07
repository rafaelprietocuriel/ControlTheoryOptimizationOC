function z = trapzsum(x,dim)
%TRAPZSUM  Trapezoidal numerical integration with unit intervals.
%
%   Z = TRAPZSUM(X,DIM) computes the integral of X along dimention DIM
%   with respect to the vector of unit intervals using the trapezoidal method.  
%     Length of vector X must be more than one ( size(X,1) > 1 ).
%
%   Z = TRAPZSUM(X) computes the integral of X along its first dimension
%   with respect to the vector of unit intervals using the trapezoidal method.  
%     Length of vector X must be more than one ( size(X,1) > 1 ).  
%
%   TRAPZSUM removes the resulting singleton.
%
%   Make sure X is an N-dimentional array with size(X,1) > 1.
%
% Anton O. Belyakov 
% Revision date: 2010/02/11

siz = size(x);
nd = length(siz); %number of dimentions (not less than 2, see help ndim)
if nargin == 2
  if dim > nd, error('Too high index Dim')
  elseif dim > 1 % put the dimention y in the first place 
    x = permute(x,[dim 1:dim-1 dim+1:nd]);
  elseif dim < 1, error('Index must be Dim > 0')
  end
end
if nd == 2   
    if siz(2) == 1 
        z = sum(x,1) - (x(1) + x(end))*0.5; % Trapezoid sum
    else
        z = (sum(x,1) - (x(1,:) + x(end,:))*0.5)'; % Trapezoid sum
    end
else
    z = sum(x,1);
    z = z(1,:) - (x(1,:) + x(end,:))*0.5; % Trapezoid sum
    z = reshape(z,siz(2:end)); % remove the first singleton dimention
end
