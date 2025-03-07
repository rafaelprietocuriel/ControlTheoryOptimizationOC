function ret=evaluatePright(y,z,coord,psival,numcols,varargin)
%Polynomial optimized version
% simplified for order one autonomous ODEs

% evaluation of the "Runge-Kutta basis polynomial" at time zero (initial
% time of the ith interval). For ODEs this value is given by y. For
% algebraic equations it is the value of the Legendre polynomial

if ~isempty(y)
    ret=y(coord,1,1);
else
    ret=sum(z(coord,:,1).*psival(ones(numel(coord),1),:,numcols+2),2);
end