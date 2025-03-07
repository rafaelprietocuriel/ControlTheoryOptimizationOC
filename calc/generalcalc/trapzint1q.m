function z = trapzint1q(x,y)
%TRAPZINT1Q  Trapezoidal numerical integration.
%
%   Z = TRAPZINT1Q(X,Y) computes the integral of Y 
%   with respect to the vector of intervals X using the trapezoidal method.  
%     Length of vector X must be one element less 
%     than the length of vector Y  ( length(X) + 1 = length(Y)).  
%
%   Make sure X and Y are column vectors.
%
% Anton O. Belyakov 
% Revision date: 2009/12/30

z = (x' * (y(1:end-1) + y(2:end)))/2;
