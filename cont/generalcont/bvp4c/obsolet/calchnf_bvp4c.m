function nfc=calchnf_bvp4c(x,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,reigvec,leigvec)

global OCMATCONT OCBVP
FcnArgs={0,freepar,modelpar};    % Pass the region index to the ODE function.

threshval=1e-6;
OCBVP.Joptions.thresh=threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh=threshval(ones(OCBVP.npar,1));


D2Phidpar2=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.nparmcod);
D2Phidpardyi=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
D2Phidyidpar=zeros(OCBVP.numode,OCBVP.numode,OCBVP.nparmcod);
D2Phidpardyip1=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
D2Phidyip1dpar=zeros(OCBVP.numode,OCBVP.numode,OCBVP.nparmcod);

D2Phidyi2=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyidyip1=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyip1dyi=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
D2Phidyip12=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);

F=OCBVP.F;



% BC points
ya=y(:,OCBVP.Lidx);
yb=y(:,OCBVP.Ridx);

rows =OCBVP.rows; 
planes=1:OCBVP.numode;

nfc=0;

%DGBC=zeros((OCBVP.numarc+1)*OCBVP.numode+OCBVP.nparmcod,1);
%dGdy=zeros(OCMATCONT.HE.numdvariablesmc,1);

% w denotes the right and v the left eigenvector in the notation of
% Kuznetsov 1998 p. 502ff
[d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=calchessianbc(bc,bcjac,bchess,ya,yb,FcnArgs{2:end});
wpar=reigvec(OCMATCONT.HE.parametermcodcoord);
reigvec=reigvec(OCMATCONT.HE.DDATA.meshvalcoord);
wa=reigvec(:,OCBVP.Lidx);
wb=reigvec(:,OCBVP.Ridx);
w1=[wa(:);wb(:);wpar(:)];
counter=0;
%row1=1:1:OCBVP.numarc*OCBVP.numode;
counterq=0;
for ii=1:OCBVP.numarc*OCBVP.numode
    counterq=counterq+1;
    counter=counter+1;
    nfc=nfc+wa(ii)*leigvec(1:OCBVP.nBCs)'*d2BCdYdya(:,:,ii)*w1;
    %DGBC(row1,1)=DGBC(row1,1)+wa(ii)*d2BCdYdya(:,:,ii)*w1;
end
counterq=0;
for ii=1:OCBVP.numarc*OCBVP.numode
    counter=counter+1;
    counterq=counterq+1;
    nfc=nfc+wb(ii)*leigvec(1:OCBVP.nBCs)'*d2BCdYdyb(:,:,ii)*w1;
    %DGBC(row1,1)=DGBC(row1,1)+wb(ii)*d2BCdYdyb(:,:,ii)*w1;
end
for ii=1:OCBVP.nparmcod
    counter=counter+1;
    nfc=nfc+w1(counter)*leigvec(1:OCBVP.nBCs)'*d2BCdYdpar(:,:,ii)*w1;
end

numarc=OCMATCONT.HE.numarc;
for arc=1:numarc
    %dGdy(planes)=dGdy(planes)+DGBC((arc-1)*OCBVP.numode+(1:OCBVP.numode));
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
    d2Fdpar2_i=d2FdYdpar_i(:,OCBVP.numode+(1:OCBVP.nparmcod),:);

    for ii=1:OCBVP.Nint(arc)
        dGdyi=zeros(OCBVP.numode,1);
        dGdyip1=zeros(OCBVP.numode,1);
        d2Fdy2_ip05mdInner_ip1=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
        d2Fdy2_ip05mdInner_i=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
        d2Fdy2_ip05mdInner_par=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
        d2Fdy2_ip05mdInner_ip1_par=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.numode);
        d2Fdy2_ip05mdInner_par_par=zeros(OCBVP.numode,OCBVP.numode,OCBVP.nparmcod);
        d2Fdy2_ip05mdInner_par_par2=zeros(OCBVP.numode,OCBVP.nparmcod,OCBVP.nparmcod);


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
                    d2Fdy2_ip05mdInner_par(:,:,kk)=d2Fdy2_ip05mdInner_par(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj+dFdy_i(ll,kk)*hi/8);
                    d2Fdy2_ip05mdInner_ip1_par(:,:,kk)=d2Fdy2_ip05mdInner_ip1_par(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj-dFdy_ip1(ll,kk)*hi/8);
                end
            end
        end
        % assembly
        for kk=1:OCBVP.nparmcod
            for ll=1:OCBVP.numode
                d2Fdy2_ip05mdInner_par_par(:,:,kk)=d2Fdy2_ip05mdInner_par_par(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
                if OCBVP.nparmcod
                    d2Fdy2_ip05mdInner_par_par2(:,:,kk)=d2Fdy2_ip05mdInner_par_par2(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
                end
            end
        end
        for jj=1:OCBVP.nparmcod
            if OCBVP.nparmcod
                D2Phidpar2(:,:,jj)=-hi/6*(d2Fdpar2_ip1(:,:,jj)-hi/2*(d2Fdy2_ip05mdInner_par_par(:,:,jj)+d2Fdydpar_ip05(:,:,jj))*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))-1/2*dFdy_ip05*(d2Fdpar2_ip1(:,:,jj)-d2Fdpar2_i(:,:,jj))*hi+4*d2Fdy2_ip05mdInner_par_par2(:,:,jj)+4*d2Fdpar2_ip05(:,:,jj)+d2Fdpar2_i(:,:,jj));
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
                D2Phidpardyi(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_i(:,:,jj)*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))*hi+1/2*dFdy_ip05*d2Fdpardy_i(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_par(:,:,jj)+d2Fdpardy_i(:,:,jj))*hi;
                D2Phidpardyip1(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_ip1(:,:,jj)*(dFdpar_ip1(:,1:OCBVP.nparmcod)-dFdpar_i(:,1:OCBVP.nparmcod))*hi-1/2*dFdy_ip05*d2Fdpardy_ip1(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_ip1_par(:,:,jj)+d2Fdpardy_ip1(:,:,jj))*hi;
            end
        end
        
        PhiHess_yi=cat(2,D2Phidyi2,D2Phidyidyip1,D2Phidpardyi);
        PhiHess_yip1=cat(2,D2Phidyip1dyi,D2Phidyip12,D2Phidpardyip1);
        PhiHessPar=cat(2,D2Phidyidpar,D2Phidyip1dpar,D2Phidpar2);
%         myy=[yi yip1];
%         numJacOpt.diffvar=2;
%         numJacOpt.vectvars=[];
%         numHess=mynumjaccsd(@jacobianPhi,{[xi xip1],myy(:),freepar,modelpar,ode,odejac},numel(myy)^2/2,numJacOpt);
%         numHess=reshape(numHess,numel(myy)/2,numel(myy),numel(myy));
%         if max(max(max(abs(numHess(:,:,1:end/2)-PhiHess_yi))))>1e-1 ||  max(max(max(abs(numHess(:,:,end/2+1:end)-PhiHess_yip1))))>1e-1
%             numHess
%         end
%         PhiHess_yi=numHess(:,:,1:end/2);
%         PhiHess_yip1=numHess(:,:,end/2+1:end);
        for jj=1:OCBVP.numode%+OCBVP.nparmcod
            counterq=counterq+1;
            nfc=nfc+wi(jj)*vt'*PhiHess_yi(:,:,jj)*wt;
            nfc=nfc+wip1(jj)*vt'*PhiHess_yip1(:,:,jj)*wt;
%             dGdyi(:,jj)=PhiHess_yi(:,:,jj)*wt;
%             dGdyip1(:,jj)=PhiHess_yip1(:,:,jj)*wt;
        end
        %dGdy(planes)=dGdy(planes)+dGdyi(1:OCBVP.numode);
        for jj=1:OCBVP.nparmcod
            nfc=nfc+wpar(jj)*vt'*PhiHessPar(:,:,jj)*wt;
        end
        %cols = cols + OCBVP.numode;
        %planes=planes+OCBVP.numode;
        %dGdy(planes)=dGdy(planes)+dGdyi(1:OCBVP.numode)+dGdyip1(1:OCBVP.numode);
        %dGdy(planes)=[dGdyi dGdyip1]*wt;

        rows = rows+OCBVP.numode;%+OCBVP.nparmcod;   % next equation
        dFdy_i=dFdy_ip1;
        dFdpar_i=dFdpar_ip1;
        d2Fdy2_i=d2Fdy2_ip1;
        d2Fdpardy_i=d2Fdpardy_ip1;
        d2Fdpar2_i=d2Fdpar2_ip1;

    end
    %dGdy(planes)=dGdy(planes)+DGBC(arc*OCBVP.numode+(1:OCBVP.numode));
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
dBCdpar=dBCdpar(:,1:OCBVP.nparmcod);

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
%dFdpar=dFdpar(:,1:OCBVP.nparmcod);
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
    reshape(D2BCdyya(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmcod,n));
dBCoptions.diffvar = 2;
D2BCdyyb=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdyb=cat(2,reshape(D2BCdyyb(1:jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyyb(jacelements+1:2*jacelements,:),OCBVP.nBCs,n,n), ...
    reshape(D2BCdyyb(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmcod,n));
dBCoptions.thresh = repmat(1e-6,npar,1);
dBCoptions.diffvar = 6;
D2BCdpar=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdpar=cat(2,reshape(D2BCdpar(1:jacelements,:),OCBVP.nBCs,n,OCBVP.npar), ...
    reshape(D2BCdpar(jacelements+1:2*jacelements,:),OCBVP.nBCs,n,OCBVP.npar), ...
    reshape(D2BCdpar(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.nparmcod,OCBVP.npar));

%---------------------------------------------------------------------------

function [dFdy,dFdp]=Fnumjac(ode,odeArgs,odeVal,Joptions,dPoptions)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.

Joptions.fac = [];
dFdy=odenumjac(ode,odeArgs,odeVal,Joptions);
dPoptions.fac = [];

dFdp=odenumjac(ode,odeArgs,odeVal,dPoptions);

%---------------------------------------------------------------------------

function [d2FdYdY,d2FdYdp]=Fnumhess(ode,odejac,x,y,n,npar,varargin)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
% Y consists of x and free parameters without the continuation parameter
% p denotes the continuation parameter

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
d2Fdpdy=d2FdYdy(end-n+1:end,:);
d2FdYdy(end-n+1:end,:)=[];
Hoptions.thresh = repmat(1e-6,npar,1);
Hoptions.diffvar = 6;
d2FdYdP=odenumjac(@odejacaux,odeArgs,odeVal,Hoptions);
d2FdYdP(end-n+1:end,:)=[];
d2FdYdp=d2FdYdP(:,end);
d2FdYdY=[d2FdYdy d2FdYdP(:,1:end-1)];
d2FdYdY=reshape(d2FdYdY,n+OCBVP.nparmcod,n+OCBVP.nparmcod,n+OCBVP.nparmcod);
d2FdYdp=reshape(d2FdYdp,n+OCBVP.nparmcod,n+OCBVP.nparmcod,1);

function J=jacobianPhi(x,y,freepar,modelpar,ode,odejac)
global OCBVP
y=reshape(y,[],2);
FcnArgs = {1,freepar,modelpar};
hreg=diff(x);
xi = x(1);
yi = y(:,1);
Fi =ode(xi,yi,FcnArgs{:});
Ji=odejac(xi,yi,FcnArgs{:});
% the right mesh point
xip1 = x(2);
yip1 = y(:,2);
Fip1 = ode(xip1,yip1,FcnArgs{:});
Jip1=odejac(xip1,yip1,FcnArgs{:});
% the midpoint
hi = hreg(1);
xip05 = (xi + xip1)/2;
yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);
Jip05 = odejac(xip05,yip05,FcnArgs{:});  % recompute the Jacobian
twiceJip05 = 2*Jip05;

J=[-(OCBVP.In+hi/6*(Ji+twiceJip05*(OCBVP.In+hi/4*Ji))) OCBVP.In-hi/6*(Jip1+twiceJip05*(OCBVP.In-hi/4*Jip1))];

function J=jacobianvect(J)

J=J(:);

function [dFdy,nfevals,nfcalls]=mynumjaccsd(fun,Fargs,nF,options)
% NUMJACCSD    Complex Step Jacobian
% based on 'jacobiancsd' by Yi Cao at Cranfield University, 02/01/2008 and
% uses the structure of odenumjac.
%
% J = NUMJACCSD(F,FARGS,NF,OPTIONS) returns the numerical (NF x N) Jacobian
% matrix of a NF-vector function, F(FARGS{:}) at the reference point, 
% Y=FARGS{OPTIONS.DIFFVAR} (N-vector). 
%
% The structure OPTIONS must have the following fields: DIFFVAR, VECTVARS.
% The field OPTIONS.DIFFVAR is the index of the  differentiation variable,
% Y = FARGS{DIFFVAR}. For a function F(t,x), set DIFFVAR to 1 to compute
% DF/Dt, or to 2 to compute DF/Dx. Set OPTIONS.VECTVAR to the indices of
% vectorized arguments: VECTVAR = [2] indicates that  F(t,[x1 y2 ...])
% returns [F(t,x1) F(t,x2) ...], while VECTVAR = [1,2] indicates that F([t1
% t2 ...],[x1 x2 ...]) returns [F(t1,x1) F(t2,x2) ...].
%   
% [DFDY,NFEVALS,NFCALLS] = ODENUMJAC(...) returns the number of values
% F(FARGS{:}) computed while forming dFdY (NFEVALS) and the number of calls
% to the function F (NFCALLS). If F is not vectorized, NFCALLS equals
% NFEVALS.  
% Dieter Grass

% Options
diffvar = options.diffvar; 
vectvar = options.vectvars;

% The differentiation variable.
y  = Fargs{diffvar};
classY = class(y);

ny = length(y);

dFdy=zeros(nF,ny);                   % allocate memory for the Jacobian matrix
del = (y + ny*eps(classY)*i) - y;
h=imag(del);
ydel = y(:,ones(1,ny)) + diag(del);
for ii=1:ny                      % loop for each independent variable
    dFdy(:,ii)=imag(jacobianvect(fun(Fargs{1:diffvar-1},ydel(:,ii),Fargs{diffvar+1:end})))/h(ii);     % complex step differentiation
end
nfcalls = ny;                       % stats

nfevals = ny;                         % stats (at least one per loop)
