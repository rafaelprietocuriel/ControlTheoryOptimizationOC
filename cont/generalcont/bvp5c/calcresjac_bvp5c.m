function Jac=calcresjac_bvp5c(X,Y,freepar,modelpar,ode,bc,odejac,bcjac)
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:54:08 $

global OCBVP OCMATCONT
ya=Y(:,1);
yb=Y(:,end);
nstages=OCBVP.nstages;
neqn=OCBVP.neqn;
neqnmc=OCBVP.neqnmc;
bigI=OCBVP.bigI;
bigA=OCBVP.bigA;
c=OCBVP.c;
threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(neqn,1));

h=diff(X(1:nstages:end));
JACbc =spalloc(neqnmc,numel(Y),2*neqn*neqn);
JACode=spalloc(numel(Y)-neqn,numel(Y),neqn*nstages*neqn*(nstages+1)*numel(h));
bcArgs={freepar,modelpar};
if isempty(bcjac)
    [dGdya,dGdyb]=BCnumjac(bc,ya,yb,bcArgs);
else
    [dGdya,dGdyb]=bcjac(ya,yb,bcArgs{:});
end
JACbc(:,1:neqn)=dGdya;
JACbc(:,end-neqn+1:end)=dGdyb;
jidx=0;
xidx=0;
Jpropagated=[];
FcnArgs={OCMATCONT.HE.arcindex,freepar,modelpar};
for i=1:numel(h)
    x=X(xidx+1:xidx+nstages+1);
    y=Y(:,xidx+1:xidx+nstages+1);
    if isempty(odejac)
        F = ode(x,y,FcnArgs{:});
        OCBVP.Fref=F(:,1);
        if isempty(Jpropagated)
            Jn=Fnumjac(ode,{x(1),y(:,1),FcnArgs{:}});
        else
            Jn=Jpropagated;
        end
        OCBVP.Fref=F(:,4);
        Jnp1=Fnumjac(ode,{x(4),y(:,4),FcnArgs{:}});
        if norm(Jnp1 - Jn,1) <= 0.25*(norm(Jnp1,1) + norm(Jn,1))
            Jnc1=(1 - c(1))*Jn + c(1)*Jnp1;
            Jnc2=(1 - c(2))*Jn + c(2)*Jnp1;
        else
            OCBVP.Fref=F(:,2);
            Jnc1=Fnumjac(ode,{x(2),y(:,2),FcnArgs{:}});
            OCBVP.Fref=F(:,3);
            Jnc2=Fnumjac(ode,{x(3),y(:,3),FcnArgs{:}});
        end
    else
        if isempty(Jpropagated)
            Jn=odejac(x(1),y(:,1),OCMATCONT.HE.arcindex,freepar,modelpar);
        else
            Jn=Jpropagated;
        end
        Jnc1=odejac(x(2),y(:,2),OCMATCONT.HE.arcindex,freepar,modelpar);
        Jnc2=odejac(x(3),y(:,3),OCMATCONT.HE.arcindex,freepar,modelpar);
        Jnp1=odejac(x(4),y(:,4),OCMATCONT.HE.arcindex,freepar,modelpar);
    end
    hJ=h(i)*[Jn,Jnc1,Jnc2,Jnp1];
    bigJ=[hJ;hJ;hJ];

    JACode( jidx+1:jidx+neqn*nstages, ...
        jidx+1:jidx+neqn*(nstages+1))=bigI + bigA.*bigJ;

    jidx=jidx + neqn*nstages;
    xidx=xidx + nstages;
end
Jac=[JACbc; JACode];