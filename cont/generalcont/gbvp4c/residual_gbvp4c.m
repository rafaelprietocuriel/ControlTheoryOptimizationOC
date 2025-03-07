function res=residual_gbvp4c(x,y,freepar,modelpar,RHS,ode)
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

%res = [];
nfcn = 0;

res=zeros(1,OCBVP.Ridx(end)-OCBVP.numarc+1);
counter=OCBVP.nBCs;
for arc = 1:OCBVP.numarc
    threshold=OCBVP.threshold(1:OCBVP.numode(arc));

    % Do not populate the interface intervals for multi-point BVPs.
    counterstep=OCBVP.numode(arc)*(OCBVP.Ridx(arc)-OCBVP.Lidx(arc));
    % Residual at the midpoints is related to the error
    % in satisfying the collocation equations.
    NewtRes=reshape(RHS(counter+(1:counterstep)),OCBVP.numode(arc),[]);
    counter=counter+counterstep;
    FcnArgs{1} = arc;

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);    % mesh point index
    Nreg = OCBVP.Nint(arc)+1;
    xreg = x(xidx);
    yreg = y(1:OCBVP.numode(arc),xidx);
    ypreg = yp(1:OCBVP.numode(arc),xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);   % mesh interval index

    thresh = threshold(:,ones(1,OCBVP.Nint(arc)));

    % the mid-points
    temp =  (NewtRes * spdiags(1.5./hreg(:),0,OCBVP.Nint(arc),OCBVP.Nint(arc))) ./ ...
        max(abs(Fmid(1:OCBVP.numode(arc),iidx)),thresh);
    res_reg = lobw3*dot(temp,temp,1);

    % Lobatto L2 points
    xLob = xreg(1:Nreg-1) + lob2*hreg;
    [yLob,ypLob] = interp_Hermite(lob2,hreg,yreg,ypreg);
    if OCMATCONT.OPTIONS.xyvectorized
        fLob = ode(xLob,yLob,FcnArgs{:});
        nfcn = nfcn + 1;
    else
        fLob = zeros(OCBVP.numode(arc),OCBVP.Nint(arc));
        for i = 1:OCBVP.Nint(arc)
            fLob(1:OCBVP.numode(arc),i) = ode(xLob(i),yLob(1:OCBVP.numode(arc),i),FcnArgs{:});
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
            fLob(1:OCBVP.numode(arc),i) = ode(xLob(i),yLob(1:OCBVP.numode(arc),i),FcnArgs{:});
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
