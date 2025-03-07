function z = trapzint3q(x,y)
%TRAPZINT3Q  Trapezoidal numerical integration.
%
%   Z = TRAPZINT3Q(X,Y) computes the integral of Y 
%   with respect to the vector of intervals X using the trapezoidal method.  
%     Length of vector X must be one element less than the length 
%     of the first dimention of array Y  ( length(X) = size(Y,1) - 1 ).  
%   TRAPZINT3Q operates along the first dimension of Y 
%     and removes the resulting singleton there.
%
%   Make sure X is a column vector and Y is an N-dimentional array (N > 2).
%
% Anton O. Belyakov 
% Revision date: 2009/12/30

% Trapezoid sum computed with vector-matrix multiply.
z = (x' * (y(1:end-1,:) + y(2:end,:)))/2;
% reshape removing the first singleton dimention
siz = size(y);
z = reshape(z,siz(2:end));
end
        

