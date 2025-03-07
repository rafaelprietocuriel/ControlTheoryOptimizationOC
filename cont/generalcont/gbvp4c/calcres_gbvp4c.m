function res=calcres_gbvp4c(x,y,freepar,modelpar,ode,bc,varargin)
%CALCRES_BVP4CV  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.
%
% adapted to the possibility that the number of ODEs is different at each
% arc and that the number of ODEs is different at each arc.

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

global OCMATCONT OCBVP

% multi-point BVP support
FcnArgs = {0,freepar,modelpar};   

F = zeros(OCBVP.maxnumode,OCBVP.N);
Fmid = zeros(OCBVP.maxnumode,OCBVP.N-1);    % include interface intervals
res = zeros(OCBVP.nBCs+sum(OCBVP.numode.*OCBVP.Nint),1);    % exclude interface intervals

% Boundary conditions
res(1:OCBVP.nBCs) = bc(y(:,OCBVP.Lidx),y(:,OCBVP.Ridx),freepar,modelpar);
%resic=0;
phiptr = OCBVP.nBCs;    % active region of res
for arc = 1:OCBVP.numarc

    FcnArgs{1} = arc;
    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);   % mesh point index
    xreg = x(xidx);
    yreg = y(1:OCBVP.numode(arc),xidx);
    Nreg=OCBVP.Nint(arc)+1;
    iidx = xidx(1:end-1);   % mesh interval index
    
    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = ode(xreg,yreg,FcnArgs{:});
    else
        Freg = zeros(OCBVP.numode(arc),Nreg);
        for i = 1:Nreg
            Freg(:,i) = ode(xreg(i),yreg(:,i),FcnArgs{:});
        end
    end

    % derivative at the midpoints
    h = diff(xreg);
    H = spdiags(h(:),0,OCBVP.Nint(arc),OCBVP.Nint(arc));
    xi = xreg(1:end-1);
    yi = yreg(:,1:end-1);
    xip1 = xreg(2:end);
    yip1 = yreg(:,2:end);
    Fi = Freg(:,1:end-1);
    Fip1 = Freg(:,2:end);

    xip05 = (xi + xip1)/2;
    yip05 = (yi + yip1)/2 - (Fip1 - Fi)*(H/8);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = ode(xip05,yip05,FcnArgs{:});
    else % not vectorized
        Fip05 = zeros(OCBVP.numode(arc),OCBVP.Nint(arc));
        for i = 1:OCBVP.Nint(arc)
            Fip05(1:OCBVP.numode(arc),i) = ode(xip05(i),yip05(:,i),FcnArgs{:});
        end
    end

    % the Lobatto IIIa formula
    Phireg = yip1 - yi - (Fip1 + 4*Fip05 + Fi)*(H/6);

    % output assembly
    res(phiptr+1:phiptr+OCBVP.numode(arc)*OCBVP.Nint(arc)) = Phireg(:);
    F(1:OCBVP.numode(arc),xidx) = Freg;
    Fmid(1:OCBVP.numode(arc),iidx) = Fip05;
    phiptr = phiptr + OCBVP.numode(arc)*OCBVP.Nint(arc);

end

OCBVP.F=F;
OCBVP.Fmid=Fmid;