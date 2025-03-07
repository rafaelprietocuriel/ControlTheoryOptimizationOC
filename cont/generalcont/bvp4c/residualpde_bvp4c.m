function res=residual_bvp4c(x,y,freepar,modelpar,RHS,ode)
global OCMATCONT OCBVP
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

Fmid=OCBVP.Fmid;
yp=OCBVP.F;
FcnArgs = {0,freepar,modelpar};   

% Lobatto quadrature
lob4 = (1+sqrt(3/7))/2;
lob2 = (1-sqrt(3/7))/2;
lobw24 = 49/90;
lobw3 = 32/45;

res = [];
nfcn = 0;
threshold=OCBVP.threshold;
% Residual at the midpoints is related to the error
% in satisfying the collocation equations.
NewtRes = zeros(OCBVP.numode,OCBVP.N-1);
% Do not populate the interface intervals for multi-point BVPs.
intidx = setdiff(1:OCBVP.N-1,OCBVP.arcposition);
NewtRes(:,intidx) = reshape(RHS(OCBVP.nBCs+OCBVP.numic+1:end),OCBVP.numode,[]);

for arc = 1:OCBVP.numarc
    FcnArgs{1} = arc;

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);    % mesh point index
    Nreg = OCBVP.Nint(arc)+1;
    xreg = x(xidx);
    yreg = y(:,xidx);
    ypreg = yp(:,xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);   % mesh interval index

    thresh = threshold(:,ones(1,OCBVP.Nint(arc)));

    % the mid-points
    temp =  (NewtRes(:,iidx) * spdiags(1.5./hreg(:),0,OCBVP.Nint(arc),OCBVP.Nint(arc))) ./ ...
        max(abs(Fmid(:,iidx)),thresh);
    res_reg = lobw3*dot(temp,temp,1);

    % Lobatto L2 points
    xLob = xreg(1:Nreg-1) + lob2*hreg;
    [yLob,ypLob] = interp_Hermite(lob2,hreg,yreg,ypreg);
    if OCMATCONT.OPTIONS.xyvectorized
        fLob = ode(xLob,yLob,FcnArgs{:});
        nfcn = nfcn + 1;
    else
        fLob = zeros(OCBVP.numode,OCBVP.Nint(arc));
        for i = 1:OCBVP.Nint(arc)
            fLob(:,i) = ode(xLob(i),yLob(:,i),FcnArgs{:});
        end
        nfcn = nfcn + OCBVP.Nint(arc);
    end
    temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
    resLob = dot(temp,temp,1);
    res_reg = res_reg + lobw24*resLob;

    % Lobatto L4 points
    xLob = xreg(1:Nreg-1) + lob4*hreg;
    [yLob,ypLob] = interp_Hermite(lob4,hreg,yreg,ypreg);
    if OCMATCONT.OPTIONS.xyvectorized
        fLob = ode(xLob,yLob,FcnArgs{:});
        nfcn = nfcn + 1;
    else
        for i = 1:OCBVP.Nint(arc)
            fLob(:,i) = ode(xLob(i),yLob(:,i),FcnArgs{:});
        end
        nfcn = nfcn + OCBVP.Nint(arc);
    end
    temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
    resLob = dot(temp,temp,1);
    res_reg = res_reg + lobw24*resLob;

    % scaling
    res_reg = sqrt( abs(hreg/2) .* res_reg);

    res(iidx) = res_reg;
end

%---------------------------------------------------------------------------

function [Sx,Spx] = interp_Hermite (hnode,diffx,y,yp)
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
