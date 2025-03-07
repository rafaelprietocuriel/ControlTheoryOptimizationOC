function res = calcres_ocbvp4c(tmesh,y,z,freepar,contval,modelpar,odefile,totalbc)
%COLLOC_RHS  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.
global OCMATCONT

% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
res=zeros(OCMATCONT.HE.numdvariables,1);
numode=OCMATCONT.DOMAINDDATA(1).numode;
numbc=sum(OCMATCONT.HE.numboundarycondition);
Fmid = zeros(numode,OCMATCONT.HE.TIMEDDATA.nummesh-1);    % include interface intervals
F = zeros(numode,OCMATCONT.HE.TIMEDDATA.nummesh);
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
        Freg = zeros(OCMATCONT.DOMAINDDATA(1).numode,Nreg);
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
    xip05 = (xi + xip1)/2;
    yip05 = (yi + yip1)/2 - (Fip1 - Fi)*(H/8);

    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = odefile(xip05,yip05,FcnArgs{:});
    else % not vectorized
        Fip05 = zeros(OCMATCONT.DOMAINDDATA(1).numode,Nint);
        for i = 1:Nint
            Fip05(:,i) = odefile(xip05(i),yip05(:,i),FcnArgs{:});
        end
    end
    % the Lobatto IIIa formula
    Phireg = yip1 - yi - (Fip1 + 4*Fip05 + Fi)*(H/6);

    % output assembly
    counterode = counterode + numode*Nint;
    res(counterode_start:counterode) = Phireg(:);

    Xmid(:,iidx) = xip05;
    Ymid(:,iidx) = yip05;
    F(:,xidx) = Freg;
    Fmid(:,iidx) = Fip05;
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