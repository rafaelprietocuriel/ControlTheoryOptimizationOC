function res=calcrespde_bvp4c(t,y,freepar,modelpar,ode,bc)
%CALCRES_BVP4C  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

global OCMATCONT OCBVP

% multi-point BVP support
FcnArgs = {0,freepar,modelpar};   

F = zeros(OCBVP.numode,OCBVP.N);
Fmid = zeros(OCBVP.numode,OCBVP.N-1);    % include interface intervals
res = zeros(OCBVP.nBCs+OCBVP.numode*(OCBVP.N-OCBVP.numarc),1);    % exclude interface intervals

% Boundary conditions
res(1:OCBVP.nBCs) = bc(y(:,OCBVP.Lidx),y(:,OCBVP.Ridx),freepar,modelpar);
if OCBVP.numic
    res(OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numic)=0;
end
resic=0;
phiptr = OCBVP.nBCs+OCBVP.numic;    % active region of res
for arc = 1:OCBVP.numarc

    FcnArgs{1} = arc;
    tidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);   % mesh point index
    treg = t(tidx);
    yreg = y(:,tidx);
    Nreg=OCBVP.Nint(arc)+1;
    iidx = tidx(1:end-1);   % mesh interval index

    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = ode(treg,yreg,FcnArgs{:});
        nfcn = 1;
    else
        Freg = zeros(OCBVP.numode,Nreg);
        for i = 1:Nreg
            Freg(:,i) = ode(treg(i),yreg(:,i),FcnArgs{:});
        end
        nfcn = Nreg;
    end

    % derivative at the midpoints
    h = diff(treg);
    H = spdiags(h(:),0,OCBVP.Nint(arc),OCBVP.Nint(arc));
    ti = treg(1:end-1);
    yi = yreg(:,1:end-1);
    tip1 = treg(2:end);
    yip1 = yreg(:,2:end);
    Fi = Freg(:,1:end-1);
    Fip1 = Freg(:,2:end);

    tip05 = (ti + tip1)/2;
    yip05 = (yi + yip1)/2 - (Fip1 - Fi)*(H/8);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = ode(tip05,yip05,FcnArgs{:});
        nfcn = nfcn + 1;
    else % not vectorized
        Fip05 = zeros(OCBVP.numode,OCBVP.Nint(arc));
        for i = 1:OCBVP.Nint(arc)
            Fip05(:,i) = ode(tip05(i),yip05(:,i),FcnArgs{:});
        end
        nfcn = nfcn + OCBVP.Nint(arc);
    end

    % the Lobatto IIIa formula
    Phireg = yip1 - yi - (Fip1 + 4*Fip05 + Fi)*(H/6);

    % output assembly
    res(phiptr+1:phiptr+OCBVP.numode*OCBVP.Nint(arc)) = Phireg(:);
    F(:,tidx) = Freg;
    Fmid(:,iidx) = Fip05;
    phiptr = phiptr + OCBVP.numode*OCBVP.Nint(arc);

end
if OCBVP.numic
    res(OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numic)=resic;
end

%OCBVP.Xmid=Xmid;
%OCBVP.Ymid=Ymid;
OCBVP.F=F;
OCBVP.Fmid=Fmid;