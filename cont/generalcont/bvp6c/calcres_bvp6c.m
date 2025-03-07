function Phi = calcres_bvp6c(x,y,freepar,modelpar,ode,bc,icfun)
%COLLOC_RHS  Evaluate the system of collocation equations Phi(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.
global OCMATCONT OCBVP

nregions = OCBVP.numarc;
Lidx = OCBVP.Lidx;
Ridx = OCBVP.Ridx;
FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.
n=OCBVP.numode;
N=OCBVP.N;

nBCs = OCBVP.nBCs;
F = zeros(n,N);
Fmid = zeros(n,N-1);    % include interface intervals
Phi = zeros(nBCs+n*(N-nregions),1);    % exclude interface intervals

% Boundary conditions
% Do not pass info about singular BVPs in ExtraArgs to BC function.
Phi(1:nBCs) = bc(y(:,Lidx),y(:,Ridx),freepar,modelpar);
if OCBVP.numic
    Phi(OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numic)=0;
end
resic=0;
phiptr = OCBVP.nBCs+OCBVP.numic;    % active region of res

for region = 1:nregions
    FcnArgs{1} = region;
    xidx = Lidx(region):Ridx(region);   % mesh point index
    Nreg = length(xidx);
    xreg = x(xidx);
    yreg = y(:,xidx);

    iidx = xidx(1:end-1);   % mesh interval index
    Nint = length(iidx);


    if OCBVP.numic
        resic=resic+icfun(xreg,yreg,FcnArgs{:});
    end
    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = ode(xreg,yreg,FcnArgs{:});
        nfcn = 1;
    else
        Freg = zeros(n,Nreg);
        for i = 1:Nreg
            Freg(:,i) = ode(xreg(i),yreg(:,i),FcnArgs{:});
        end
        nfcn = Nreg;
    end

    %mesh point data
    h = diff(xreg);
    H = spdiags(h(:),0,Nint,Nint);
    xi = xreg(1:end-1);
    yi = yreg(:,1:end-1);
    xip1 = xreg(2:end);
    yip1 = yreg(:,2:end);
    Fi = Freg(:,1:end-1);
    Fip1 = Freg(:,2:end);

    %interior points & derivative
    xip025 = 0.25*(3*xi + xip1);
    xip075 = 0.25*(xi + 3*xip1);
    yip025 = (54*yi + 10*yip1 + (9*Fi - 3*Fip1)*H)/64;
    yip075 = (10*yi + 54*yip1 + (3*Fi - 9*Fip1)*H)/64;
    if OCMATCONT.OPTIONS.xyvectorized
        Fip025 = ode(xip025,yip025,FcnArgs{:});
        Fip075 = ode(xip075,yip075,FcnArgs{:});
        nfcn = nfcn + 3;
    else % not vectorized
        Fip025 = zeros(n,Nint);
        Fip075 = zeros(n,Nint);
        for i = 1:Nint
            Fip025(:,i) = ode(xip025(i),yip025(:,i),FcnArgs{:});
            Fip075(:,i) = ode(xip075(i),yip075(:,i),FcnArgs{:});
        end
        nfcn = nfcn + 2*Nint;
    end

    %mid points & derivative
    xip05 = 0.5*(xi + xip1);
    yip05 = 0.5*(yi + yip1) - (5*Fi - 16*Fip025 + 16*Fip075 - 5*Fip1)*(H/24);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = ode(xip05,yip05,FcnArgs{:});
        nfcn = nfcn + 1;
    else % not vectorized
        Fip05 = zeros(n,Nint);
        for i = 1:Nint
            Fip05(:,i) = ode(xip05(i),yip05(:,i),FcnArgs{:});
        end
        nfcn = nfcn + Nint;
    end

    % the Cash-Singhal formula
    Phireg = yip1 - yi - (7*Fi + 32*Fip025 + 12*Fip05 + 32*Fip075 + 7*Fip1)*(H/90);

    % output assembly
    Phi(phiptr+1:phiptr+n*Nint) = Phireg(:);
    phiptr = phiptr + n*Nint;

    Xmid(:,iidx,1) = xip025;
    Xmid(:,iidx,2) = xip05;
    Xmid(:,iidx,3) = xip075;

    Ymid(:,iidx,1) = yip025;
    Ymid(:,iidx,2) = yip05;
    Ymid(:,iidx,3) = yip075;

    F(:,xidx) = Freg;
    Fmid(:,iidx,1) = Fip025;
    Fmid(:,iidx,2) = Fip05;
    Fmid(:,iidx,3) = Fip075;

end
if OCBVP.numic
    Phi(OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numic)=resic;
end

OCBVP.Xmid=Xmid;
OCBVP.Ymid=Ymid;
OCBVP.F=F;
OCBVP.Fmid=Fmid;
