function res = calcres_ocbvp6c(tmesh,y,z,freepar,contval,modelpar,odefile,bc)
%COLLOC_RHS  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.
global OCMATCONT

% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
counterode=sum(OCMATCONT.HE.numboundarycondition);
res=zeros(OCMATCONT.HE.numdvariables,1);

numode=OCMATCONT.DOMAINDDATA(1).numode;
F = zeros(numode,OCMATCONT.HE.TIMEDDATA.nummesh);
Fmid = zeros(numode,OCMATCONT.HE.TIMEDDATA.nummesh-1,3);    % include interface intervals
% Boundary conditions
FcnArgs={modelpar,freepar,contval};
res(1:numbc)=totalbc(y(:,leftarcindex),y(:,rightarcindex),FcnArgs{:});
% Do not pass info about singular BVPs in ExtraArgs to BC function.
counterode=numbc;
for arc=1:numarc
    counterode_start=counterode+1;%numode+OCMATCONT.HE.numparametermc;
    FcnArgs={modelpar,freepar,contval,arc};
    xidx = leftarcindex(arc):rightarcindex(arc);   % mesh point index
    Nreg = length(xidx);
    xreg = tmesh(xidx);
    yreg = y(:,xidx);

    iidx = xidx(1:end-1);   % mesh interval index
    Nint = length(iidx);


    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = odefile(xreg,yreg,FcnArgs{:});
    else
        Freg = zeros(numode,Nreg);
        for i = 1:Nreg
            Freg(:,i) = odefile(xreg(i),yreg(:,i),FcnArgs{:});
        end
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
        Fip025 = odefile(xip025,yip025,FcnArgs{:});
        Fip075 = odefile(xip075,yip075,FcnArgs{:});
    else % not vectorized
        Fip025 = zeros(numode,Nint);
        Fip075 = zeros(numode,Nint);
        for i = 1:Nint
            Fip025(:,i) = odefile(xip025(i),yip025(:,i),FcnArgs{:});
            Fip075(:,i) = odefile(xip075(i),yip075(:,i),FcnArgs{:});
        end
    end
    %mid points & derivative
    xip05 = 0.5*(xi + xip1);
    yip05 = 0.5*(yi + yip1) - (5*Fi - 16*Fip025 + 16*Fip075 - 5*Fip1)*(H/24);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = odefile(xip05,yip05,FcnArgs{:});
    else % not vectorized
        Fip05 = zeros(numode,Nint);
        for i = 1:Nint
            Fip05(:,i) = odefile(xip05(i),yip05(:,i),FcnArgs{:});
        end
    end
    % the Cash-Singhal formula
    Phireg = yip1 - yi - (7*Fi + 32*Fip025 + 12*Fip05 + 32*Fip075 + 7*Fip1)*(H/90);

    % output assembly
    counterode = counterode + numode*Nint;
    res(counterode_start:counterode) = Phireg(:);

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

if OCMATCONT.storedata
    OCMATCONT.HE.DDATA.Xmid=Xmid;
    OCMATCONT.HE.DDATA.Ymid=Ymid;
    OCMATCONT.HE.DDATA.F=F;
    OCMATCONT.HE.DDATA.Fmid=Fmid;
end

function [y,unknownPar]=reshapeY(Y,n,N,nN,npar)

y=reshape(Y(1:nN),n,N);
unknownPar{1}=Y(nN+1:nN+npar);