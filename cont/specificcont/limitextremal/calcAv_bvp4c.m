function dGdy=calcAv_bvp4c(x,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,reigvec)

global OCMATCONT OCBVP
FcnArgs={0,freepar,modelpar};    % Pass the region index to the ODE function.

threshval=1e-6;
OCBVP.Joptions.thresh=threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh=threshval(ones(OCBVP.npar,1));

D2Phidpardyi=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
D2Phidyidpar=zeros(OCBVP.numode,OCBVP.numode,OCBVP.npar);
D2Phidpardyip1=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
D2Phidyip1dpar=zeros(OCBVP.numode,OCBVP.numode,OCBVP.npar);

D2Phidyi2=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyidyip1=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyip1dyi=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyip12=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidpar2=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.npar);
F=OCBVP.F;

wpar=reigvec(OCMATCONT.HE.parametermcodcoord);
reigvec=reigvec(OCMATCONT.HE.DDATA.meshvalcoord);
% BC points
ya=y(:,OCBVP.Lidx);
yb=y(:,OCBVP.Ridx);

rows =OCBVP.rows;   % define the action area
cols =1:OCBVP.numode;             % in the global Jacobian
planes=1:OCBVP.numode;
par_planesidx=OCMATCONT.HE.parametercoord;
dGdy=zeros(OCMATCONT.HE.numdvariables,1);
%par_rowidx=OCMATCONT.HE.numdvariables-OCMATCONT.HE.numparameter+(1:OCMATCONT.HE.numparameter-1);

[d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=calchessianbc(bc,bcjac,bchess,ya,yb,FcnArgs{2:end});
wa=reigvec(:,OCBVP.Lidx);
wb=reigvec(:,OCBVP.Ridx);
w1=[wa(:);wb(:);wpar(:)];
counter=0;
DGBC=zeros((OCBVP.numarc+1)*OCBVP.numode+OCBVP.nparmcod,1);
dGtotpar=zeros(OCBVP.npar,1);

for ii=1:OCBVP.numarc*OCBVP.numode
    counter=counter+1;
    DGBC(counter,1)=-v(1:OCBVP.nBCs)'*d2BCdYdya(:,:,ii)*w1;
end
for ii=1:OCBVP.numarc*OCBVP.numode
    counter=counter+1;
    DGBC(counter,1)=-v(1:OCBVP.nBCs)'*d2BCdYdyb(:,:,ii)*w1;
end
for ii=1:OCBVP.npar
    dGtotpar(ii,1)=-v(1:OCBVP.nBCs)'*d2BCdYdpar(:,:,ii)*w1;
end
dGdy(par_planesidx)=dGtotpar;
% DGBC is the derivative of the equation G for the minimally extended
% system, to continue a limit point solution, with respect to ya, vb and
% the free parameters p 


numarc=OCMATCONT.HE.numarc;
dGdyi=zeros(OCBVP.numode,1);
dGdyip1=zeros(OCBVP.numode,1);
for arc=1:numarc
    dGdy(planes)=dGdy(planes)+DGBC((arc-1)*OCBVP.numode+(1:OCBVP.numode));
    FcnArgs{1}=arc;

    xidx=OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    xreg=x(xidx);
    yreg=y(:,xidx);
    wreg=reigvec(:,xidx);
    Freg=F(:,xidx);
    hreg=diff(xreg);
    

    [dFdy_i dFdpar_i] =calcjacobianode(ode,odejac,xreg(1),yreg(:,1),FcnArgs{:});
    [d2FdYdy_i d2FdYdpar_i]=calchessianode(ode,odejac,odehess,xreg(1),yreg(:,1),FcnArgs{:});
    d2Fdy2_i=d2FdYdy_i(:,1:OCBVP.numode,:);
    d2Fdpardy_i=d2FdYdy_i(:,OCBVP.numode+(1:OCBVP.nparmcod),:);
    d2Fdydpar_i=d2FdYdpar_i(:,1:OCBVP.numode,:);
    d2Fdpar2_i=d2FdYdpar_i(:,OCBVP.numode+(1:OCBVP.nparmcod),:);
    %nfcn=nfcn+nFcalls;
    for ii=1:OCBVP.Nint(arc)

        d2Fdy2_ip05mdInner_ip1=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
        d2Fdy2_ip05mdInner_i=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
        d2Fdy2_ip05mdInner_par=zeros(OCBVP.numode,OCBVP.npar,OCBVP.numode);
        d2Fdy2_ip05mdInner_ip1_par=zeros(OCBVP.numode,OCBVP.npar,OCBVP.numode);
        d2Fdy2_ip05mdInner_par_par=zeros(OCBVP.numode,OCBVP.numode,OCBVP.npar);
        d2Fdy2_ip05mdInner_par_par2=zeros(OCBVP.numode,OCBVP.npar,OCBVP.npar);


        % the left mesh point
        xi=xreg(ii);
        yi=yreg(:,ii);
        wi=wreg(:,ii);
        Fi=Freg(:,ii);
        % the right mesh point
        xip1=xreg(ii+1);
        yip1=yreg(:,ii+1);
        wip1=wreg(:,ii+1);
        Fip1=Freg(:,ii+1);
        OCBVP.Fref=Fip1;
        [dFdy_ip1 dFdpar_ip1]=calcjacobianode(ode,odejac,xip1,yip1,FcnArgs{:});
        [d2FdYdy_ip1 d2FdYdpar_ip1]=calchessianode(ode,odejac,odehess,xip1,yip1,FcnArgs{:});
        hi=hreg(ii);
        xip05=(xi + xip1)/2;
        yip05=(yi+yip1)/2-hi/8*(Fip1-Fi);
        [dFdy_ip05 dFdpar_ip05]=calcjacobianode(ode,odejac,xip05,yip05,FcnArgs{:});
        [d2FdYdy_ip05 d2FdYdpar_ip05]=calchessianode(ode,odejac,odehess,xip05,yip05,FcnArgs{:});

        d2Fdy2_ip1=d2FdYdy_ip1(:,1:OCBVP.numode,:);
        d2Fdpardy_ip1=d2FdYdy_ip1(:,OCBVP.numode+(1:OCBVP.nparmcod),:);
        d2Fdydpar_ip1=d2FdYdpar_ip1(:,1:OCBVP.numode,:);
        d2Fdpar2_ip1=d2FdYdpar_ip1(:,OCBVP.numode+(1:OCBVP.nparmcod),:);

        % the midpoint
        d2Fdy2_ip05=d2FdYdy_ip05(:,1:OCBVP.numode,:);
        d2Fdpardy_ip05=d2FdYdy_ip05(:,OCBVP.numode+(1:OCBVP.nparmcod),:);
        d2Fdydpar_ip05=d2FdYdpar_ip05(:,1:OCBVP.numode,:);
        d2Fdpar2_ip05=d2FdYdpar_ip05(:,OCBVP.numode+(1:OCBVP.nparmcod),:);

        wt=[wi;wip1;wpar];
        vt=leigvec(rows);
        for kk=1:OCBVP.numode
            for ll=1:OCBVP.numode
                Iniijj=OCBVP.In(kk,ll)/2;
                d2Fdy2_ip05mdInner_ip1(:,:,kk)=d2Fdy2_ip05mdInner_ip1(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(Iniijj-dFdy_ip1(ll,kk)*hi/8);
                d2Fdy2_ip05mdInner_i(:,:,kk)=d2Fdy2_ip05mdInner_i(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(Iniijj+dFdy_i(ll,kk)*hi/8);
                if OCBVP.nparmcod
                    d2Fdy2_ip05mdInner_par(:,1:OCBVP.nparmcod,kk)=d2Fdy2_ip05mdInner_par(:,1:OCBVP.nparmcod,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj+dFdy_i(ll,kk)*hi/8);
                    d2Fdy2_ip05mdInner_ip1_par(:,1:OCBVP.nparmcod,kk)=d2Fdy2_ip05mdInner_ip1_par(:,1:OCBVP.nparmcod,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj-dFdy_ip1(ll,kk)*hi/8);
                end
            end
        end
        % assembly
        for kk=1:OCBVP.npar
            for ll=1:OCBVP.numode
                d2Fdy2_ip05mdInner_par_par(:,:,kk)=d2Fdy2_ip05mdInner_par_par(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
                if OCBVP.nparmcod
                    d2Fdy2_ip05mdInner_par_par2(:,1:OCBVP.nparmcod,kk)=d2Fdy2_ip05mdInner_par_par2(:,1:OCBVP.nparmcod,kk)+d2Fdpardy_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
                end
            end
        end
        for jj=1:OCBVP.npar
            if OCBVP.nparmcod
                D2Phidpar2(:,:,jj)=-hi/6*(d2Fdpar2_ip1(:,:,jj)-hi/2*(d2Fdy2_ip05mdInner_par_par(:,:,jj)+d2Fdydpar_ip05(:,:,jj))*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))-1/2*dFdy_ip05*(d2Fdpar2_ip1(:,:,jj)-d2Fdpar2_i(:,:,jj))*hi+4*d2Fdy2_ip05mdInner_par_par2(:,1:OCBVP.nparmcod,jj)+4*d2Fdpar2_ip05(:,:,jj)+d2Fdpar2_i(:,:,jj));
            end
            D2Phidyidpar(:,:,jj)=-1/6*(4*((d2Fdy2_ip05mdInner_par_par(:,:,jj)+d2Fdydpar_ip05(:,:,jj))*(OCBVP.In/2+1/8*dFdy_i*hi)+1/8*dFdy_ip05*d2Fdydpar_i(:,:,jj)*hi)+d2Fdydpar_i(:,:,jj))*hi;
            D2Phidyip1dpar(:,:,jj)=-1/6*(4*((d2Fdy2_ip05mdInner_par_par(:,:,jj)+d2Fdydpar_ip05(:,:,jj))*(OCBVP.In/2-1/8*dFdy_ip1*hi)-1/8*dFdy_ip05*d2Fdydpar_ip1(:,:,jj)*hi)+d2Fdydpar_ip1(:,:,jj))*hi;
        end

        for jj=1:OCBVP.numode
            D2Phidyi2(:,:,jj)=-1/6*(2*(d2Fdy2_ip05mdInner_i(:,:,jj)*(OCBVP.In+1/4*dFdy_i*hi)+1/4*dFdy_ip05*d2Fdy2_i(:,:,jj)*hi)+d2Fdy2_i(:,:,jj))*hi;

            D2Phidyidyip1(:,:,jj)=-1/3*d2Fdy2_ip05mdInner_i(:,:,jj)*(OCBVP.In-1/4*dFdy_ip1*hi)*hi;
            D2Phidyip1dyi(:,jj,1:OCBVP.numode)=D2Phidyidyip1(:,:,jj);
            
            D2Phidyip12(:,:,jj)=-1/6*(d2Fdy2_ip1(:,:,jj)+2*(d2Fdy2_ip05mdInner_ip1(:,:,jj)*(OCBVP.In-1/4*dFdy_ip1*hi)-1/4*dFdy_ip05*d2Fdy2_ip1(:,:,jj)*hi))*hi;
            if OCBVP.nparmcod
                D2Phidpardyi(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_i(:,:,jj)*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))*hi+1/2*dFdy_ip05*d2Fdpardy_i(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_par(:,1:OCBVP.nparmcod,jj)+d2Fdpardy_i(:,:,jj))*hi;
                D2Phidpardyip1(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_ip1(:,:,jj)*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))*hi-1/2*dFdy_ip05*d2Fdpardy_ip1(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_ip1_par(:,1:OCBVP.nparmcod,jj)+d2Fdpardy_ip1(:,:,jj))*hi;
            end
            %D2Phidyidpar(:,jj,1:OCBVP.npar)=D2Phidpardyi(:,:,jj);
            %D2Phidyip1dpar(:,jj,1:OCBVP.npar)=D2Phidpardyip1(:,:,jj);
        end
        PhiHess_yi=cat(2,D2Phidyi2,D2Phidyidyip1,D2Phidpardyi);
        PhiHess_yip1=cat(2,D2Phidyip1dyi,D2Phidyip12,D2Phidpardyip1);
        PhiHessPar=cat(2,D2Phidyidpar,D2Phidyip1dpar,D2Phidpar2);

        for jj=1:OCBVP.numode
            dGdyi(jj)=-vt'*PhiHess_yi(:,:,jj)*wt;
            dGdyip1(jj)=-vt'*PhiHess_yip1(:,:,jj)*wt;
        end
        dGdy(planes)=dGdy(planes)+dGdyi(1:OCBVP.numode);
        for jj=1:OCBVP.npar
            dGtotpar(jj)=-vt'*PhiHessPar(:,:,jj)*wt;
        end
        dGdy(par_planesidx)=dGdy(par_planesidx)+dGtotpar;

        cols = cols + OCBVP.numode;
        planes=planes+OCBVP.numode;
        dGdy(planes)=dGdy(planes)+dGdyip1(1:OCBVP.numode);

        %last_planes(rows,cols,:)=D2Phidyidpar;
        rows = rows+OCBVP.numode;   % next equation
        dFdy_i=dFdy_ip1;
        dFdpar_i=dFdpar_ip1;
        d2Fdy2_i=d2Fdy2_ip1;
        d2Fdpar2_i=d2Fdpar2_ip1;
    end
    dGdy(planes)=dGdy(planes)+DGBC(arc*OCBVP.numode+(1:OCBVP.numode));
end

function [d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=calchessianbc(bc,bcjac,bchess,ya,yb,varargin)
global OCBVP
if isempty(bchess)
    [d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=BCnumhess(bc,bcjac,ya,yb,OCBVP.numode,OCBVP.npar,varargin{:});
else
    [d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=bchess(ya,yb,varargin{:});
end

function [dBCdya,dBCddyb,dBCdpar]=calcjacobianbc(bc,bcjac,ya,yb,varargin)
global OCBVP

if isempty(bcjac)
    [dBCdya,dBCddyb,dBCdpar]=BCnumjac(bc,ya,yb,OCBVP.numode,OCBVP.npar,varargin{:});
else
    [dBCdya,dBCddyb,dBCdpar]=bcjac(ya,yb,varargin{:});
end
dBCdpar=dBCdpar(:,1:OCBVP.nparmc);

function [d2FdYdy,d2FdYdpar]=calchessianode(ode,odejac,odehess,x,y,varargin)
global OCBVP

if isempty(odehess)
    [d2FdYdy,d2FdYdpar]=Fnumhess(ode,odejac,x,y,OCBVP.numode,OCBVP.npar,varargin{:});
else
    [d2FdYdy,d2FdYdpar]=odehess(x,y,varargin{:});
end

function [dFdy,dFdpar]=calcjacobianode(ode,odejac,x,y,varargin)
global OCBVP

OCBVP.Joptions.vectvars=[];
if isempty(odejac)
    odeVal=ode(x,y,varargin{:});
    [dFdy,dFdpar]=Fnumjac(ode,{x,y,varargin{:}},odeVal,OCBVP.Joptions,OCBVP.dPoptions);
else
    [dFdy,dFdpar]=odejac(x,y,varargin{:});
end
%dFdpar=dFdpar(:,1:OCBVP.nparmc);
%---------------------------------------------------------------------------

function res = bcaux(Ya,Yb,n,bcfun,varargin)
ya = reshape(Ya,n,[]);
yb = reshape(Yb,n,[]);
res = bcfun(ya,yb,varargin{:});

%---------------------------------------------------------------------------

function resjac = bcjacaux(Ya,Yb,n,bc,bcjac,varargin)
ya = reshape(Ya,n,[]);
yb = reshape(Yb,n,[]);
[dBCdya,dBCddyb,dBCdpar]=calcjacobianbc(bc,bcjac,ya,yb,varargin{:});
resjac=[dBCdya(:);dBCddyb(:);dBCdpar(:)];

%---------------------------------------------------------------------------

function jac = odejacaux(x,y,ode,odejac,varargin)
% dFdpar: par minus continuation parameter
[dFdy dFdpar]=calcjacobianode(ode,odejac,x,y,varargin{:});
jac=[dFdy(:);dFdpar(:)];

%---------------------------------------------------------------------------

function [dBCdya,dBCdyb,dBCdpar] = BCnumjac(bc,ya,yb,n,npar,varargin)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.

% Do not pass info about singular BVPs in ExtraArgs to BC function.
bcArgs = {ya(:),yb(:),n,bc,varargin{:}};
dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
dBCoptions.fac = [];
dBCoptions.vectvars = []; % BC functions not vectorized

bcVal  = bcaux(bcArgs{:});

dBCoptions.diffvar = 1;
dBCdya=odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
dBCoptions.diffvar = 2;
dBCdyb=odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
bcArgs = {ya,yb,varargin{:}};
dBCoptions.thresh = repmat(1e-6,npar,1);
dBCoptions.diffvar = 3;
dBCdpar=odenumjac(bc,bcArgs,bcVal,dBCoptions);

%---------------------------------------------------------------------------

function [d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=BCnumhess(bc,bcjac,ya,yb,n,npar,varargin)
global OCBVP

% Do not pass info about singular BVPs in ExtraArgs to BC function.
bcArgs = {ya(:),yb(:),n,bc,bcjac,varargin{:}};
dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
dBCoptions.fac = [];
dBCoptions.vectvars = []; % BC functions not vectorized

bcVal  = bcjacaux(bcArgs{:});
jacelements=n*OCBVP.nBCs;
dBCoptions.diffvar = 1;
D2BCdyya=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdya=cat(2,reshape(D2BCdyya(1:jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyya(jacelements+1:2*jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyya(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmc,n));
dBCoptions.diffvar = 2;
D2BCdyyb=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdyb=cat(2,reshape(D2BCdyyb(1:jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyyb(jacelements+1:2*jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyyb(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmc,n));
dBCoptions.thresh = repmat(1e-6,npar,1);
dBCoptions.diffvar = 6;
D2BCdpar=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdpar=cat(2,reshape(D2BCdpar(1:jacelements,:),OCBVP.nBCs,n,OCBVP.npar), ...
    reshape(D2BCdpar(jacelements+1:2*jacelements,:),OCBVP.nBCs,n,OCBVP.npar), ...
    reshape(D2BCdpar(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmc,OCBVP.npar));

%---------------------------------------------------------------------------

function [dFdy,dFdp]=Fnumjac(ode,odeArgs,odeVal,Joptions,dPoptions)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.

Joptions.fac = [];
dFdy=odenumjac(ode,odeArgs,odeVal,Joptions);
dPoptions.fac = [];

dFdp=odenumjac(ode,odeArgs,odeVal,dPoptions);

%---------------------------------------------------------------------------

function [d2FdYdy,d2FdYdp]=Fnumhess(ode,odejac,x,y,n,npar,varargin)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
global OCBVP
if isempty(odejac)
    warning('Results may not be accurate, due to numerical difficulties in the calculation of the Hessian')
end
Hoptions.thresh = repmat(1e-6,n,1);
Hoptions.fac = [];
Hoptions.vectvars = []; % BC functions not vectorized
Hoptions.diffvar = 2;
odeArgs = {x,y,ode,odejac,varargin{:}};

odeVal=odejacaux(x,y,ode,odejac,varargin{:});
d2FdYdy=odenumjac(@odejacaux,odeArgs,odeVal,Hoptions);
d2FdYdy=reshape(d2FdYdy,n+OCBVP.npar,n,n);
Hoptions.thresh = repmat(1e-6,npar,1);
Hoptions.diffvar = 6;
d2FdYdp=odenumjac(@odejacaux,odeArgs,odeVal,Hoptions);
d2FdYdp=reshape(d2FdYdp,n+OCBVP.npar,n,npar);
