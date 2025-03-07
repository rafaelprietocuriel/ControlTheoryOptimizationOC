function [vHw vHwpar]=calchessianphi_ocbvp4c(tmesh,y,z,freepar,contval,modelpar,equationfile,jacobianfile,jacobianfilepar,bc,bcjacobianfile,hessianfile,hessianfilepar,bchessianfile,w,v)

global OCMATCONT OCMATLSC
domainddata=OCMATCONT.DOMAINDDATA;

if 1%nargin==14
    w=reshape(OCMATLSC.LAE_phi,size(y,1),[]);
    v=reshape(OCMATLSC.LAE_psi,size(y,1),[]);
else
    w=reshape(w,size(y,1),[]);
    v=reshape(v,size(y,1),[]);
end
F=OCMATCONT.HE.DDATA.F;
% multi-point BVP support
numarc=OCMATCONT.HE.numarc;

leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
In=eye(domainddata.numode);
% numeric calculation of boundary Hessian
[numrow numcolumn]=size(OCMATCONT.HE.yL);
nphase=2*numrow*numcolumn;
yL=OCMATCONT.HE.yL(:);
yR=OCMATCONT.HE.yR(:);
FcnArgs={modelpar,freepar,contval};
HBC=zeros(OCMATCONT.HE.totalnumboundarycondition,2*numcolumn*numrow+2,2*numcolumn*numrow+2);
ybc=[yL(:);yR(:)];
for ii=1:nphase
    ybc1=ybc;
    ybc1(ii)=ybc1(ii)-OCMATCONT.OPTIONS.increment;
    ybc2=ybc;
    ybc2(ii)=ybc2(ii)+OCMATCONT.OPTIONS.increment;
    ybc2=reshape(ybc2,numrow,2*numcolumn);
    ybc1=reshape(ybc1,numrow,2*numcolumn);
    HBC(:,:,ii)=jacbc(ybc2(:,1:numcolumn),ybc2(:,numcolumn+1:2*numcolumn),bcjacobianfile,FcnArgs)-jacbc(ybc1(:,1:numcolumn),ybc1(:,numcolumn+1:2*numcolumn),bcjacobianfile,FcnArgs);
end
ybc=reshape(ybc,numrow,2*numcolumn);
for ii=1:2
    contval1=contval;
    contval1(ii)=contval1(ii)-OCMATCONT.OPTIONS.increment;
    contval2=contval;
    contval1(ii)=contval1(ii)+OCMATCONT.OPTIONS.increment;
    FcnArgs1={modelpar,freepar,contval1};
    FcnArgs2={modelpar,freepar,contval2};
    HBC(:,:,nphase+ii)=jacbc(ybc(:,1:numcolumn),ybc(:,numcolumn+1:2*numcolumn),bcjacobianfile,FcnArgs2)-jacbc(ybc1(:,1:numcolumn),ybc1(:,numcolumn+1:2*numcolumn),bcjacobianfile,FcnArgs1);
end
HBC=HBC/(2*OCMATCONT.OPTIONS.increment);
HBC(:,nphase+(1:2),:)=[];
wL=w(:,leftarcindex);
wR=w(:,rightarcindex);
vL=v(:,leftarcindex);
HBCpar=HBC(:,:,nphase+(1:2));
vHw=[vL(:,1)']*[tprod(HBC(:,:,1:nphase/2)+HBC(:,:,nphase/2+1:nphase),[1 -1 2],[wL(:);wR(:)],-1)];
vHwpar=[vL(:,1)']*[tprod(HBCpar,[1 -1 2],[wL(:);wR(:)],-1)];
vHwa=vL(:,1)'*tprod(HBC(:,:,1:nphase/2),[1 -1 2],[wL(:);wR(:)],-1);
vHwb=vL(:,1)'*tprod(HBC(:,:,nphase/2+1:nphase),[1 -1 2],[wL(:);wR(:)],-1);
vHw=zeros(1,numel(y));
vHw(1:numrow)=vHwa;
counter=0;
counter_start=counter+1;
counter=counter+nphase/2;
for arc=1:numarc
    FcnArgs={modelpar,freepar,contval,arc};
    FcnArgspar={modelpar,freepar,arc};

    xidx = leftarcindex(arc):rightarcindex(arc);
    xreg = tmesh(xidx);
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);
    wreg=w(:,xidx);
    vreg=v(:,xidx);
    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);

    % Collocation equations
    Ji=calcjacobian(equationfile,jacobianfile,xreg(1),yreg(:,1),FcnArgs,1);
    Ji_par=calcjacobianp(equationfile,jacobianfilepar,xreg(1),yreg(:,1),contval,FcnArgspar,1);
    Hi=calchessian(equationfile,jacobianfile,hessianfile,xreg(1),yreg(:,1),FcnArgs,1);
    Hi_par=calchessianp(equationfile,jacobianfile,hessianfilepar,xreg(1),yreg(:,1),contval,FcnArgspar,1);
    for ii = 1:Nint
        hi = hreg(ii);
        % the left mesh point
        xi = xreg(ii);
        yi = yreg(:,ii);
        wi = wreg(:,ii);
        Fi = Freg(:,ii);
        % the right mesh point
        xip1 = xreg(ii+1);
        yip1 = yreg(:,ii+1);
        wip1=wreg(:,ii+1);
        vip1=vreg(:,ii+1);
        Fip1 = Freg(:,ii+1);
        Jip1=calcjacobian(equationfile,jacobianfile,xip1,yip1,FcnArgs,1);
        Hip1=calchessian(equationfile,jacobianfile,hessianfile,xip1,yip1,FcnArgs,1);
        Jip1_par=calcjacobianp(equationfile,jacobianfilepar,xip1,yip1,contval,FcnArgspar,1);
        Hip1_par=calchessianp(equationfile,jacobianfile,hessianfilepar,xip1,yip1,contval,FcnArgspar,1);

        %the interior points
        xip05 = (xi + xip1)/2;
        yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);

        Hip05=calchessian(equationfile,jacobianfile,hessianfile,xip05,yip05,FcnArgs,1);
        Jip05=calcjacobian(equationfile,jacobianfile,xip05,yip05,FcnArgs,1);
        Jip05_par=calcjacobianp(equationfile,jacobianfilepar,xip05,yip05,contval,FcnArgspar,1);
        Hip05_par=calchessianp(equationfile,jacobianfile,hessianfilepar,xip05,yip05,contval,FcnArgspar,1);

       
        D2PHIY02=hi*(-1/6*Hi-2/12*tprod(tprod(Hip05,[1 2 -1],(In+1/4*hi*Ji),[-1 3]),[1 -1 3],(In+1/4*hi*Ji),[-1 2])-1/12*hi*tprod(Jip05,[1 -1],Hi,[-1 2 3]));
        D2PHIY0Y1=-2/3*hi*tprod(tprod(Hip05,[1 2 -1],(1/2*In-1/8*hi*Jip1),[-1 3]),[1 -1 3],(1/2*In+1/8*hi*Ji),[-1 2]);
        D2PHIY1Y0=-2/3*hi*tprod(tprod(Hip05,[1 2 -1],(1/2*In+1/8*hi*Ji),[-1 3]),[1 -1 3],(1/2*In-1/8*hi*Jip1),[-1 2]);
        D2PHIY12=hi*(-1/6*Hip1-2/12*tprod(tprod(Hip05,[1 2 -1],(In-1/4*hi*Jip1),[-1 3]),[1 -1 3],(In-1/4*hi*Jip1),[-1 2])+1/12*hi*tprod(Jip05,[1 -1],Hip1,[-1 2 3]));

        D2PHIDY0DP=hi*(-1/6*Hi_par-2/3*tprod((hi/8*tprod(Hip05,[1 2 -1],(Ji_par-Jip1_par),[-1 3])+Hip05_par),[1 -1 3],(1/2*In+1/8*hi*Ji),[-1 2])-hi/12*tprod(Jip05,[1 -1],Hi_par,[-1 2 3]));
        D2PHIDY1DP=hi*(-1/6*Hip1_par-2/3*tprod((hi/8*tprod(Hip05,[1 2 -1],(Ji_par-Jip1_par),[-1 3])+Hip05_par),[1 -1 3],(1/2*In-1/8*hi*Jip1),[-1 2])-hi/12*tprod(Jip05,[-1 1],Hip1_par,[-1 2 3]));

%         yt=[yi;yip1];
%         nphase=numel(yt);
%         for jj=1:nphase
%             yt1=yt;
%             yt1(jj)=yt1(jj)-OCMATCONT.OPTIONS.increment;
%             yt2=yt;
%             yt2(jj)=yt2(jj)+OCMATCONT.OPTIONS.increment;
%             D2PHIY2(:,:,jj)=jacphi([xi xip1],reshape(yt2,nphase/2,2),equationfile,jacobianfile,hi,In,FcnArgs)-jacphi([xi xip1],reshape(yt1,nphase/2,2),equationfile,jacobianfile,hi,In,FcnArgs);
%         end
%         D2PHIY2=D2PHIY2/(2*OCMATCONT.OPTIONS.increment);

        %vHw=[vHw vip1'*[tprod(D2PHIY02,[1 -1 2],wi,-1)+tprod(D2PHIY1Y0,[1 -1 2],wip1,-1)+tprod(D2PHIY0Y1,[1 -1 2],wi,-1)+tprod(D2PHIY12,[1 -1 2],wip1,-1)]];
        vHwi=vip1'*[tprod(D2PHIY02,[1 -1 2],wi,-1)+tprod(D2PHIY1Y0,[1 -1 2],wip1,-1)];
        vHw(counter_start:counter)=vHw(counter_start:counter)+vHwi;
        vHwip1=vip1'*[tprod(D2PHIY0Y1,[1 -1 2],wi,-1)+tprod(D2PHIY12,[1 -1 2],wip1,-1)];
        counter_start=counter+1;
        counter=counter+nphase/2;
        vHw(counter_start:counter)=vHwip1;
        vHwpar=vHwpar+vip1'*[tprod(D2PHIDY0DP,[1 -1 2],wi,-1)+tprod(D2PHIDY1DP,[1 -1 2],wip1,-1)];
        Ji = Jip1;
        Hi = Hip1;
    end
    vHw(counter_start:counter)=vHw(counter_start:counter)+vHwb;
end


function h=calchessian(equationfile,jacobianfile,hessianfile,t,x,FcnArgs,modelflag)
%
% h=chess(equationfile,jacobianfile,hessianfile,x,p)
% Calculates the hessian of the system
%
global OCMATCONT

if modelflag
    if OCMATCONT.OPTIONS.SymDerivative >= 2
        % Use symbolic derivatives if they are defined
        h = feval(hessianfile,t,x,FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(x,1);
        for ii=1:nphase
            x1 = x; x1(ii) = x1(ii)-OCMATCONT.OPTIONS.increment;
            x2 = x; x2(ii) = x2(ii)+OCMATCONT.OPTIONS.increment;
            h(:,:,ii)=calcjacobian(equationfile,jacobianfile,t,x2,FcnArgs,modelflag)-calcjacobian(equationfile,jacobianfile,t,x1,FcnArgs,modelflag);
        end
        h = h/(2*OCMATCONT.OPTIONS.increment);
    end

else
    if OCMATCONT.symhess
        % Use symbolic derivatives if they are defined
        h = feval(hessianfile,x);
    else
        % If not, use finite differences
        for ii=1:OCMATCONT.ndim
            x1 = x; x1(ii) = x1(ii)-OCMATCONT.OPTIONS.increment;
            x2 = x; x2(ii) = x2(ii)+OCMATCONT.OPTIONS.increment;
            h(:,:,ii)=calcjacobian(equationfile,jacobianfile,x2,[])-calcjacobian(equationfile,jacobianfile,x1,[]);
        end
        h = h/(2*OCMATCONT.OPTIONS.increment);
    end
end

function h=calchessianp(equationfile,jacobianfile,hessianfile,t,x,p,FcnArgs,modelflag)
%
% h=chess(equationfile,jacobianfile,hessianfile,x,p)
% Calculates the hessian of the system
%
global OCMATCONT
FcnArgsnew([1 2 4])=FcnArgs;
if modelflag
    if OCMATCONT.OPTIONS.SymDerivative >= 2
        % Use symbolic derivatives if they are defined
        h = feval(hessianfile,t,x,p,FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(p,1);
        for ii=1:nphase
            p1 = p; p1(ii) = p1(ii)-OCMATCONT.OPTIONS.increment;
            p2 = p; p2(ii) = p2(ii)+OCMATCONT.OPTIONS.increment;
            FcnArgs1=FcnArgsnew;
            FcnArgs1{3}=p1;
            FcnArgs2=FcnArgsnew;
            FcnArgs2{3}=p2;
            h(:,:,ii)=calcjacobian(equationfile,jacobianfile,t,x,FcnArgs1,modelflag)-calcjacobian(equationfile,jacobianfile,t,x,FcnArgs2,modelflag);
        end
        h = h/(2*OCMATCONT.OPTIONS.increment);
    end

else
    if OCMATCONT.symhess
        % Use symbolic derivatives if they are defined
        h = feval(hessianfile,x,p);
    else
        % If not, use finite differences
        for ii=1:OCMATCONT.ndim
            x1 = x; x1(ii) = x1(ii)-OCMATCONT.OPTIONS.increment;
            x2 = x; x2(ii) = x2(ii)+OCMATCONT.OPTIONS.increment;
            h(:,:,ii)=calcjacobianp(equationfile,jacobianfile,x2,[])-calcjacobianp(equationfile,jacobianfile,x1,[]);
        end
        h = h/(2*OCMATCONT.OPTIONS.increment);
    end
end

function j=calcjacobian(equationfile,jacobianfile,t,x,FcnArgs,modelflag)
%
% j=cjac(equationfile,jacobianfile,x,p)
% Calculates the jacobianfile of the system.
% p can be empty, if the active parameter is included in x (and the correct
% Jacobian function is passed).
%
global OCMATCONT

if modelflag==1
    if (OCMATCONT.OPTIONS.SymDerivative >=1)
        % Use symbolic derivatives if they are defined
        j = feval(jacobianfile,t,x, FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(x,1);
        for ii=1: nphase
            x1 = x; x1(ii) = x1(ii)-OCMATCONT.OPTIONS.increment;
            x2 = x; x2(ii) = x2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(t,x2,FcnArgs{:})-equationfile(t,x1,FcnArgs{:});
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
elseif modelflag==2 % without time argument
    if (OCMATCONT.OPTIONS.SymDerivative >=1)
        % Use symbolic derivatives if they are defined
        j = feval(jacobianfile,x,FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(x,1);
        for ii=1: nphase
            x1 = x; x1(ii) = x1(ii)-OCMATCONT.OPTIONS.increment;
            x2 = x; x2(ii) = x2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(x2,FcnArgs{:})-equationfile(x1,FcnArgs{:});
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
else
    if OCMATCONT.symjac
        % Use curve-derivatives if they are defined
        j = feval(jacobianfile,x);
    else
        % If not, use finite differences
        x1 = x;
        x2 = x;
        for ii=1:OCMATCONT.ndim
            x1(ii)=x1(ii)-OCMATCONT.OPTIONS.increment;
            x2(ii)=x2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(t,x2)-equationfile(t,x1);
            x1(ii)=x(ii);
            x2(ii)=x(ii);
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
end

function j=calcjacobianp(equationfile,jacobianfile,t,x,p,FcnArgs,modelflag)
%
% j=cjac(equationfile,jacobianfile,x,p)
% Calculates the jacobianfile of the system.
% p can be empty, if the active parameter is included in x (and the correct
% Jacobian function is passed).
%
global OCMATCONT

if modelflag==1
    if (OCMATCONT.OPTIONS.SymDerivative >=1)
        % Use symbolic derivatives if they are defined
        j = feval(jacobianfile,t,x,p,FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(p,1);
        for ii=1: nphase
            p1 = p; p1(ii) = p1(ii)-OCMATCONT.OPTIONS.increment;
            p2 = p; p2(ii) = p2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(t,x,p,FcnArgs{:})-equationfile(t,x,p,FcnArgs{:});
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
elseif modelflag==2 % without time argument
    if (OCMATCONT.OPTIONS.SymDerivative >=1)
        % Use symbolic derivatives if they are defined
        j = feval(jacobianfile,x,p,FcnArgs{:});
    else
        % If not, use finite differences
        nphase=size(p,1);
        for ii=1: nphase
            p1 = p; p1(ii) = p1(ii)-OCMATCONT.OPTIONS.increment;
            p2 = p; p2(ii) = p2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(x,p,FcnArgs{:})-equationfile(x,p,FcnArgs{:});
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
else
    if OCMATCONT.symjac
        % Use curve-derivatives if they are defined
        j = feval(jacobianfile,x,p);
    else
        % If not, use finite differences
        x1 = x;
        x2 = x;
        for ii=1:OCMATCONT.ndim
            x1(ii)=x1(ii)-OCMATCONT.OPTIONS.increment;
            x2(ii)=x2(ii)+OCMATCONT.OPTIONS.increment;
            j(:,ii)=equationfile(t,x2)-equationfile(t,x1);
            x1(ii)=x(ii);
            x2(ii)=x(ii);
        end
        j = j/(2*OCMATCONT.OPTIONS.increment);
    end
end

function [J]=jacbc(yL,yR,bcjacobianfile,FcnArgs)
[Ja Jb Jpar]=bcjacobianfile(yL,yR,FcnArgs{:});

J=[Ja Jb Jpar];
