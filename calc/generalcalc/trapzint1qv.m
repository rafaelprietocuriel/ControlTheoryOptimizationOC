function z = trapzint1qv(x,y)
%TRAPZINT1Q  Trapezoidal numerical integration.
%
%   Z = TRAPZINT1QV(X,Y) computes the integral of Y 
%   with respect to the vector of intervals X using the trapezoidal method.  
%     Length of vector X must be one element less 
%     than the length of vector Y  ( length(X) + 1 = length(Y)).  
%
%   Make sure X and Y are column vectors.
%
% Anton O. Belyakov 
% Revision date: 2009/12/30
z = sum(x(:,ones(1,size(y,3))).* (y(1:end-1,:) + y(2:end,:))).'/2;
