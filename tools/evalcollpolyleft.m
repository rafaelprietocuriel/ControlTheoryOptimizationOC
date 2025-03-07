function ret=evalcollpolyleft(ii,y,z,coord,psival,numcols,aeflag)
%Polynomial optimized version
% simplified for order one autonomous ODEs

% evaluation of the "Runge-Kutta basis polynomial" at time zero (initial
% time of the ith interval). For ODEs this value is given by y. For
% algebraic equations it is the value of the Legendre polynomial

if ~aeflag
    ret=y(coord,1,ii);
else
    ret=sum(z(coord,:,ii).*psival(ones(numel(coord),1),:,numcols+2),2);
end