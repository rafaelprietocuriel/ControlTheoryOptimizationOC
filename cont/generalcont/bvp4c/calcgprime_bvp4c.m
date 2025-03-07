function dGdy=calcgprime_bvp4c(x,y,freepar,modelpar,ode,bc,odejac,bcjac,odehess,bchess,leigvec,reigvec)

global OCMATCONT OCBVP
F=OCBVP.F;

FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.

% BC points
ya = y(:,OCBVP.Lidx);
yb = y(:,OCBVP.Ridx);


%%%%
%% HESSIAN
%%%%
zerosInnerTemplate=zeros(OCBVP.numode,OCBVP.numode,OCBVP.numode);
zerosInnerTemplate2=zeros(OCBVP.numode,OCBVP.npar,OCBVP.numode);
zerosInnerTemplate3=zeros(OCBVP.numode,OCBVP.numode,OCBVP.npar);
zerosInnerTemplate4=zeros(OCBVP.numode,OCBVP.npar,OCBVP.npar);

D2Phidpardyi=zerosInnerTemplate2;
D2Phidyidpar=zerosInnerTemplate3;
D2Phidpardyip1=zerosInnerTemplate2;
D2Phidyip1dpar=zerosInnerTemplate3;

D2Phidyi2=zerosInnerTemplate;
D2Phidyidyip1=zerosInnerTemplate;
D2Phidyip1dyi=zerosInnerTemplate;
D2Phidyip12=zerosInnerTemplate;
D2Phidpar2=zerosInnerTemplate4;

planes=1:OCBVP.numode;
rows = OCBVP.rows;   % define the action area
cols = OCBVP.cols;             % in the global Jacobian
planespar=OCBVP.nN+(1:OCBVP.npar);
planesparmc=OCBVP.nN+(1:OCBVP.nparmc);
colsparmc=OCBVP.nN+(1:OCBVP.nparmc);

dGdy=zeros(OCBVP.nN+OCBVP.npar,1);
[d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=calchessianbc(bc,bcjac,bchess,ya,yb,FcnArgs{2:end});

d2BCdyadya=d2BCdYdya(:,1:OCBVP.numode*OCMATCONT.HE.numarc,:);
d2BCdybdya=d2BCdYdya(:,OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.numode*OCMATCONT.HE.numarc),:);
d2BCdpardya=d2BCdYdya(:,2*OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.npar),:);
d2BCdyadyb=d2BCdYdyb(:,1:OCBVP.numode*OCMATCONT.HE.numarc,:);
d2BCdybdyb=d2BCdYdyb(:,OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.numode*OCMATCONT.HE.numarc),:);
d2BCdpardyb=d2BCdYdyb(:,2*OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.npar),:);
d2BCdyadpar=d2BCdYdpar(:,1:OCBVP.numode*OCMATCONT.HE.numarc,:);
d2BCdybdpar=d2BCdYdpar(:,OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.numode*OCMATCONT.HE.numarc),:);
d2BCdpardpar=d2BCdYdpar(:,2*OCBVP.numode*OCMATCONT.HE.numarc+(1:OCBVP.npar),:);

numarc=OCMATCONT.HE.numarc;
for arc=1:numarc
    colsa=cols;
    planesa=planes;
    FcnArgs{1}=arc;

    xidx=OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    xreg=x(xidx);
    yreg=y(:,xidx);
    Freg=F(:,xidx);
    hreg=diff(xreg);
    

    [dFdy_i,dFdpar_i]=calcjacobianode(ode,odejac,xreg(1),yreg(:,1),FcnArgs{:});
    [d2Fdy2_i,d2Fdydpar_i,d2Fdpar2_i]=calchessianode(ode,odejac,odehess,xreg(1),yreg(:,1),FcnArgs{:});
    counter=0;
    for ii=1:OCBVP.Nint(arc)
        counter=counter+1;

        d2Fdy2_ip05mdInner_ip1=zerosInnerTemplate;
        d2Fdy2_ip05mdInner_i=zerosInnerTemplate;
        d2Fdy2_ip05mdInner_i_par=zerosInnerTemplate2;
        d2Fdy2_ip05mdInner_ip1_par=zerosInnerTemplate2;
        d2Fdy2_ip05mdInner_par_par=zerosInnerTemplate3;
        d2Fdy2_ip05mdInner_par_par2=zerosInnerTemplate4;

        % the left mesh point
        xi=xreg(ii);
        yi=yreg(:,ii);
        Fi=Freg(:,ii);
        % the right mesh point
        xip1=xreg(ii+1);
        yip1=yreg(:,ii+1);
        Fip1=Freg(:,ii+1);
        OCBVP.Fref=Fip1;
        [dFdy_ip1,dFdpar_ip1]=calcjacobianode(ode,odejac,xip1,yip1,FcnArgs{:});
        [d2Fdy2_ip1,d2Fdydpar_ip1 d2Fdpar2_ip1]=calchessianode(ode,odejac,odehess,xip1,yip1,FcnArgs{:});
        hi=hreg(ii);
        xip05=(xi + xip1)/2;
        yip05=(yi+yip1)/2-hi/8*(Fip1-Fi);
        [dFdy_ip05,dFdpar_ip05]=calcjacobianode(ode,odejac,xip05,yip05,FcnArgs{:});
        [d2Fdy2_ip05,d2Fdydpar_ip05,d2Fdpar2_ip05]=calchessianode(ode,odejac,odehess,xip05,yip05,FcnArgs{:});
        d2Fdpardy_ip05=permute(d2Fdydpar_ip05,[1 3 2]);
        d2Fdpardy_i=permute(d2Fdydpar_i,[1 3 2]);
        d2Fdpardy_ip1=permute(d2Fdydpar_ip1,[1 3 2]);

        for kk=1:OCBVP.numode
            for ll=1:OCBVP.numode
                Iniijj=OCBVP.In(kk,ll)/2;
                d2Fdy2_ip05mdInner_i(:,:,kk)=d2Fdy2_ip05mdInner_i(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(Iniijj+dFdy_i(ll,kk)*hi/8);
                d2Fdy2_ip05mdInner_ip1(:,:,kk)=d2Fdy2_ip05mdInner_ip1(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(Iniijj-dFdy_ip1(ll,kk)*hi/8);
                
                d2Fdy2_ip05mdInner_i_par(:,:,kk)=d2Fdy2_ip05mdInner_i_par(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj+dFdy_i(ll,kk)*hi/8);
                d2Fdy2_ip05mdInner_ip1_par(:,:,kk)=d2Fdy2_ip05mdInner_ip1_par(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(Iniijj-dFdy_ip1(ll,kk)*hi/8);
            end
        end
        % assembly
        for kk=1:OCBVP.npar
            for ll=1:OCBVP.numode
                d2Fdy2_ip05mdInner_par_par(:,:,kk)=d2Fdy2_ip05mdInner_par_par(:,:,kk)+d2Fdy2_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
                d2Fdy2_ip05mdInner_par_par2(:,:,kk)=d2Fdy2_ip05mdInner_par_par2(:,:,kk)+d2Fdpardy_ip05(:,:,ll)*(dFdpar_i(ll,kk)-dFdpar_ip1(ll,kk))*hi/8;
            end
        end
        for jj=1:OCBVP.numode
            D2Phidyi2(:,:,jj)=-1/6*(2*(d2Fdy2_ip05mdInner_i(:,:,jj)*(OCBVP.In+1/4*dFdy_i*hi)+1/4*dFdy_ip05*d2Fdy2_i(:,:,jj)*hi)+d2Fdy2_i(:,:,jj))*hi;

            D2Phidyip1dyi(:,:,jj)=-1/3*d2Fdy2_ip05mdInner_i(:,:,jj)*(OCBVP.In-1/4*dFdy_ip1*hi)*hi;
            D2Phidyidyip1(:,jj,:)=D2Phidyip1dyi(:,:,jj);

            D2Phidyip12(:,:,jj)=-1/6*(d2Fdy2_ip1(:,:,jj)+2*(d2Fdy2_ip05mdInner_ip1(:,:,jj)*(OCBVP.In-1/4*dFdy_ip1*hi)-1/4*dFdy_ip05*d2Fdy2_ip1(:,:,jj)*hi))*hi;
            
            D2Phidpardyi(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_i(:,:,jj)*(dFdpar_ip1(:,:)-dFdpar_i(:,:))*hi+1/2*dFdy_ip05*d2Fdpardy_i(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_i_par(:,:,jj)+d2Fdpardy_i(:,:,jj))*hi;
            D2Phidpardyip1(:,:,jj)=-1/6*(-1/2*d2Fdy2_ip05mdInner_ip1(:,:,jj)*(dFdpar_ip1(:,:)-dFdpar_i(:,:))*hi-1/2*dFdy_ip05*d2Fdpardy_ip1(:,:,jj)*hi+4*d2Fdy2_ip05mdInner_ip1_par(:,:,jj)+d2Fdpardy_ip1(:,:,jj))*hi;
            
            D2Phidyidpar(:,jj,:)=D2Phidpardyi(:,:,jj);
            D2Phidyip1dpar(:,jj,:)=D2Phidpardyip1(:,:,jj);
        end
        for jj=1:OCBVP.npar
            D2Phidpar2(:,:,jj)=-hi/6*(d2Fdpar2_ip1(:,:,jj)-hi/2*(d2Fdy2_ip05mdInner_par_par(:,:,jj)+d2Fdydpar_ip05(:,:,jj))*(dFdpar_ip1(:,:)-dFdpar_i(:,:))-1/2*dFdy_ip05*(d2Fdpar2_ip1(:,:,jj)-d2Fdpar2_i(:,:,jj))*hi+4*d2Fdy2_ip05mdInner_par_par2(:,:,jj)+4*d2Fdpar2_ip05(:,:,jj)+d2Fdpar2_i(:,:,jj));
        end
        PhiHess_yi=cat(2,D2Phidyi2,D2Phidyip1dyi);
        PhiHess_yip1=cat(2,D2Phidyidyip1,D2Phidyip12);
        PhiHessPar=cat(2,D2Phidyidpar,D2Phidyip1dpar);
        reigvec_c=reigvec([cols cols + OCBVP.numode]);
        leigvec_r=leigvec(rows);
        for pp=1:OCBVP.numode
            dGdy(planes(pp))=dGdy(planes(pp))-leigvec_r'*PhiHess_yi(:,:,pp)*reigvec_c;
        end

        %Hess(rows,[cols cols + OCBVP.numode],planes)=PhiHess_yi;
        planes=planes+OCBVP.numode;
        for pp=1:OCBVP.numode
            dGdy(planes(pp))=dGdy(planes(pp))-leigvec_r'*PhiHess_yip1(:,:,pp)*reigvec_c;
        end

        %Hess(rows,[cols cols + OCBVP.numode],planes)=PhiHess_yip1;
        if OCBVP.npar>0
            for pp=1:OCBVP.npar
                dGdy(planespar(pp))=dGdy(planespar(pp))-leigvec_r'*PhiHessPar(:,:,pp)*reigvec_c;
            end
            PhiHessParT=permute(PhiHessPar,[1 3 2]);
            reigvec_c=reigvec(planesparmc);
            pidx=[cols cols + OCBVP.numode];
            for pp=1:2*OCBVP.numode
                dGdy(pidx(pp))=dGdy(pidx(pp))-leigvec_r'*PhiHessParT(:,1:OCBVP.nparmc,pp)*reigvec_c;
            end
            reigvec_c=reigvec(planesparmc);
            for pp=1:OCBVP.npar
                dGdy(planespar(pp))=dGdy(planespar(pp))-leigvec_r'*D2Phidpar2(:,1:OCBVP.nparmc,pp)*reigvec_c;
            end
        end
        %last_planes(rows,[cols cols + OCBVP.numode],:)=PhiHessPar;
        %last_planes2(rows,1:OCBVP.npar,:)=D2Phidpar2;
        cols = cols + OCBVP.numode;

        rows = rows+OCBVP.numode;   % next equation
        dFdy_i=dFdy_ip1;
        dFdpar_i=dFdpar_ip1;
        
        d2Fdy2_i=d2Fdy2_ip1;
        d2Fdydpar_i=d2Fdydpar_ip1;
        d2Fdpar2_i=d2Fdpar2_ip1;
    end
    colsb=cols;
    planesb=planes;
    arcidx=(arc-1)*OCBVP.numode+(1:OCBVP.numode);
    leigvec_r=leigvec(1:OCBVP.nBCs);
    if OCBVP.nparmc>0
        reigvec_c=reigvec(colsa);
        for pp=1:OCBVP.numode
            dGdy(planesa(pp))=dGdy(planesa(pp))-leigvec_r'*d2BCdyadya(:,arcidx,arcidx(pp))*reigvec_c;
        end
        
        reigvec_c=reigvec(colsb);
        for pp=1:OCBVP.numode
            dGdy(planesa(pp))=dGdy(planesa(pp))-leigvec_r'*d2BCdybdya(:,arcidx,arcidx(pp))*reigvec_c;
        end
        
        reigvec_c=reigvec(colsparmc);
        for pp=1:OCBVP.numode
            dGdy(planesa(pp))=dGdy(planesa(pp))-leigvec_r'*d2BCdpardya(:,1:OCBVP.nparmc,arcidx(pp))*reigvec_c;
        end
        
        reigvec_c=reigvec(colsa);
        for pp=1:OCBVP.numode
            dGdy(planesb(pp))=dGdy(planesb(pp))-leigvec_r'*d2BCdyadyb(:,arcidx,arcidx(pp))*reigvec_c;
        end
        
        reigvec_c=reigvec(colsb);
        for pp=1:OCBVP.numode
            dGdy(planesb(pp))=dGdy(planesb(pp))-leigvec_r'*d2BCdybdyb(:,arcidx,arcidx(pp))*reigvec_c;
        end
        
        reigvec_c=reigvec(colsparmc);
        for pp=1:OCBVP.numode
            dGdy(planesb(pp))=dGdy(planesb(pp))-leigvec_r'*d2BCdpardyb(:,1:OCBVP.nparmc,arcidx(pp))*reigvec_c;
        end

        
        reigvec_c=reigvec(colsa);
        for pp=1:OCBVP.nparmc
            dGdy(planesparmc(pp))=dGdy(planesparmc(pp))-leigvec_r'*d2BCdyadpar(:,arcidx,pp)*reigvec_c;
        end
        
        reigvec_c=reigvec(colsb);
        for pp=1:OCBVP.nparmc
            dGdy(planesparmc(pp))=dGdy(planesparmc(pp))-leigvec_r'*d2BCdybdpar(:,arcidx,pp)*reigvec_c;
        end
        
        reigvec_c=reigvec(colsparmc);
        for pp=1:OCBVP.nparmc
            dGdy(planesparmc(pp))=dGdy(planesparmc(pp))-leigvec_r'*d2BCdpardpar(:,1:OCBVP.nparmc,pp)*reigvec_c;
        end
        
    end
    %     Hess(1:OCBVP.nBCs,colsa,planesa) = d2BCdyadya(:,arcidx,arcidx);
    %     Hess(1:OCBVP.nBCs,colsb,planesa) = d2BCdybdya(:,arcidx,arcidx);
    %     Hess(1:OCBVP.nBCs,colspar,planesa) = d2BCdpardya(:,1:OCBVP.npar,arcidx);
    %     Hess(1:OCBVP.nBCs,colsa,planesb) = d2BCdyadyb(:,arcidx,arcidx);
    %     Hess(1:OCBVP.nBCs,colsb,planesb) = d2BCdybdyb(:,arcidx,arcidx);
    %     Hess(1:OCBVP.nBCs,colspar,planesb) = d2BCdpardyb(:,1:OCBVP.npar,arcidx);
    %     Hess(1:OCBVP.nBCs,colsa,planespar) = d2BCdyadpar(:,arcidx,1:OCBVP.npar);
    %     Hess(1:OCBVP.nBCs,colsb,planespar) = d2BCdybdpar(:,arcidx,1:OCBVP.npar);
    %     Hess(1:OCBVP.nBCs,colspar,planespar) = d2BCdpardpar(:,1:OCBVP.npar,1:OCBVP.npar);
    cols = cols + OCBVP.numode;
    planes=planes+OCBVP.numode;
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
    [dBCdya,dBCddyb,dum,dBCdpar]=BCnumjac(bc,ya,yb,varargin);
else
    [dBCdya,dBCddyb,dBCdpar]=bcjac(ya,yb,varargin{:});
end
%dBCdpar=dBCdpar(:,1:OCBVP.nparmc);

function [d2FdY2,d2FdYdpar,d2Fdpar2]=calchessianode(ode,odejac,odehess,x,y,varargin)
global OCBVP

if isempty(odehess)
    [d2FdY2,d2FdYdpar,d2Fdpar2]=Fnumhess(ode,odejac,x,y,OCBVP.numode,OCBVP.npar,varargin{:});
else
    [d2FdY2,d2FdYdpar,d2Fdpar2]=odehess(x,y,varargin{:});
    %[d2FdY2test,d2FdYdpartest,d2Fdpar2test]=Fnumhess(ode,odejac,x,y,OCBVP.numode,OCBVP.npar,varargin{:});
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

%---------------------------------------------------------------------------

function resjac = bcjacaux(Ya,Yb,n,bc,bcjac,varargin)
ya = reshape(Ya,n,[]);
yb = reshape(Yb,n,[]);
[dBCdya,dBCddyb,dBCdpar]=calcjacobianbc(bc,bcjac,ya,yb,varargin{:});
resjac=[dBCdya(:);dBCddyb(:);dBCdpar(:)];

%---------------------------------------------------------------------------

function jac = odejacaux(x,y,ode,odejac,varargin)
[dFdy dFdpar]=calcjacobianode(ode,odejac,x,y,varargin{:});
jac=[dFdy(:);dFdpar(:)];

%---------------------------------------------------------------------------

function [d2BCdYdya,d2BCdYdyb,d2BCdYdpar]=BCnumhess(bc,bcjac,ya,yb,n,npar,varargin)
global OCBVP

% Do not pass info about singular BVPs in ExtraArgs to BC function.
bcArgs = {ya(:),yb(:),n,bc,bcjac,varargin{:}};
dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
dBCoptions.fac = [];
dBCoptions.vectvars = []; % BC functions not vectorized

freepar=varargin{1};
modelpar=varargin{2};
bcVal  = bcjacaux(bcArgs{:});
jacelements=OCBVP.numarc*n*OCBVP.nBCs;
dBCoptions.diffvar = 1;

D2BCdyya=numjaccsd(@bcjacaux,bcArgs,OCBVP.nBCs*(2*OCBVP.numarc*n+length(freepar)),dBCoptions);
%D2BCdyyaTest=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdya=cat(2,reshape(D2BCdyya(1:jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.numarc*n), ...
    reshape(D2BCdyya(jacelements+1:2*jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.numarc*n), ...
    reshape(D2BCdyya(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.npar,OCBVP.numarc*n));
dBCoptions.diffvar = 2;
D2BCdyyb=numjaccsd(@bcjacaux,bcArgs,OCBVP.nBCs*(2*OCBVP.numarc*n+length(freepar)),dBCoptions);
%D2BCdyyb=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdyb=cat(2,reshape(D2BCdyyb(1:jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.numarc*n), ...
    reshape(D2BCdyyb(jacelements+1:2*jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.numarc*n), ...
    reshape(D2BCdyyb(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.npar,OCBVP.numarc*n));
dBCoptions.thresh = repmat(1e-6,npar,1);
dBCoptions.diffvar = 6;
D2BCdpar=numjaccsd(@bcjacaux,bcArgs,OCBVP.nBCs*(2*OCBVP.numarc*n+length(freepar)),dBCoptions);
%D2BCdpar=odenumjac(@bcjacaux,bcArgs,bcVal,dBCoptions);
d2BCdYdpar=cat(2,reshape(D2BCdpar(1:jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.npar), ...
    reshape(D2BCdpar(jacelements+1:2*jacelements,:),OCBVP.nBCs,OCBVP.numarc*n,OCBVP.npar), ...
    reshape(D2BCdpar(2*jacelements+1:end,:),OCBVP.nBCs,OCBVP.npar,OCBVP.npar));

%---------------------------------------------------------------------------

function [dFdy,dFdp]=Fnumjac(ode,odeArgs,odeVal,Joptions,dPoptions)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.

Joptions.fac = [];
dFdy=odenumjac(ode,odeArgs,odeVal,Joptions);
dPoptions.fac = [];

dFdp=odenumjac(ode,odeArgs,odeVal,dPoptions);

%---------------------------------------------------------------------------

function [d2FdY2,d2FdYdp,d2Fdp2]=Fnumhess(ode,odejac,x,y,n,npar,varargin)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
global OCBVP

Hoptions.vectvars = []; % BC functions not vectorized
Hoptions.diffvar = 2;
odeArgs = {x,y,ode,odejac,varargin{:}};
d2FdYdy=numjaccsd(@odejacaux,odeArgs,n*(n+OCBVP.npar),Hoptions);
d2FdY2=reshape(d2FdYdy(1:n^2,:),n,n,n);
%d2FdY2=d2FdYdy(1:n,1:n,1:n);
Hoptions.diffvar = 6;
d2FdYdP=numjaccsd(@odejacaux,odeArgs,n*(n+OCBVP.npar),Hoptions);
d2FdYdp=reshape(d2FdYdP(1:n^2,:),n,n,npar);
d2Fdp2=reshape(d2FdYdP(n^2+1:end,:),n,npar,npar);


function Hessnum=numhess(x,y,freepar,modelpar)

options.diffvar=2;
options.vectvars=[];
Hessnum=numjaccsd(@jacaux,{x,y(:),freepar,modelpar},(numel(y)+length(freepar)-1)*(numel(y)+length(freepar)),options);
Hessnum=reshape(Hessnum,numel(y)+length(freepar)-1,numel(y)+length(freepar),[]);
options.diffvar=3;
options.vectvars=[];
Hessnumpar=numjaccsd(@jacaux,{x,y(:),freepar,modelpar},(numel(y)+length(freepar)-1)*(numel(y)+length(freepar)),options);
Hessnum=cat(3,Hessnum,reshape(Hessnumpar,numel(y)+length(freepar)-1,numel(y)+length(freepar),[]));
%---------------------------------------------------------------------------

function Jac = jacaux(x,y,freepar,modelpar)
% dFdpar: par minus continuation parameter
y=reshape(y,[],length(x));
Jac=DPhiglobal(x,y,freepar,modelpar);
Jac=Jac(:);