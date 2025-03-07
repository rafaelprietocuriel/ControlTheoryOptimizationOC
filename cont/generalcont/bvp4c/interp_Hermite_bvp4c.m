function [Sx,Spx] = interp_Hermite_bvp4c(hnode,diffx,y,yp)
%INTERP_HERMITE  Evaluate the cubic Hermite interpolant and its first 
%   derivative at x+hnode*diffx. 
N = size(y,2);
diffx = diffx(:);  % column vector
diagscal = spdiags(1./diffx,0,N-1,N-1);
slope = (y(:,2:N) - y(:,1:N-1)) * diagscal;
c = (3*slope - 2*yp(:,1:N-1) - yp(:,2:N)) * diagscal;
d = (yp(:,1:N-1)+yp(:,2:N)-2*slope) * (diagscal.^2);

diagscal = spdiags(hnode*diffx,0,diagscal);
d = d*diagscal;

Sx = ((d+c)*diagscal+yp(:,1:N-1))*diagscal+y(:,1:N-1);
Spx = (3*d+2*c)*diagscal+yp(:,1:N-1);
