function ret=evaluatePleft(ii,diffmesh,y,z,coord,psival,numcols,varargin)
% simplified for order one autonomous ODEs

% evaluation of the "Runge-Kutta basis polynomial" at time one (end time of
% the ith interval). For ODEs this value is given by y. For algebraic
% equations it is the value of the Legendre polynomial.

if ~isempty(y)
    ret=y(coord,1,ii)+sum(diffmesh(ii)*z(coord,:,ii).*psival(ones(numel(coord),1),:,numcols+1),2);
else
    ret=sum(z(coord,:,ii).*psival(ones(numel(coord),1),:,numcols+1),2);
end