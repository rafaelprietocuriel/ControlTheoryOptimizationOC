function z = trapzint2q(x,y)
%TRAPZINT2Q  Trapezoidal numerical integration.
%
%   Z = TRAPZINT2Q(X,Y) computes the integral of Y 
%   with respect to the vector of intervals X using the trapezoidal method.  
%     Length of vector X must be one element less 
%     than the row number of matrix Y  ( length(X) = size(Y,1) - 1 ).  
%   TRAPZINT2Q operates along the first dimension of Y 
%     and removes the resulting singleton there by trasposition.
%
%   Make sure X and Y are column vectors, or Y is a matrix.
%
% Anton O. Belyakov 
% Revision date: 2009/12/30

% Trapezoid sum computed with vector-matrix multiply.
z = ((y(1:end-1,:) + y(2:end,:))' * x)/2;
