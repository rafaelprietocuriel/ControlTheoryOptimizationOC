function [J,ROW,COL]=calcresdynamicjac_bvp6c(tmesh,y,z,freepar,contval,modelpar,jacobianfile,rowoffset)

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

Xmid=OCMATCONT.HE.DDATA.Xmid;
Ymid=OCMATCONT.HE.DDATA.Ymid;
F=OCMATCONT.HE.DDATA.F;
% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;

nonzeroent=OCMATCONT.HE.numdvariables*(2*domainddata(1).numeq);
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
J=zeros(nonzeroent,1);
In=eye(domainddata.numode);

gen_c=0; %General count for sparse matrix efficiency
ROWcounter=rowoffset;
for arc=1:numarc
    FcnArgs={modelpar,freepar,contval,arc};

    arcindex=OCMATCONT.HE.arcindex(arc);
    numode=domainddata(arcindex).numode;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
    
    COLcounter=0;
    xidx = leftarcindex(arc):rightarcindex(arc);
    xreg = tmesh(xidx);
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);

    [X1qtrreg, Xmidreg, X3qtrreg] = midptreg(iidx,Xmid);
    [Y1qtrreg, Ymidreg, Y3qtrreg] = midptreg(iidx,Ymid);
    %[F1qtrreg, Fmidreg, F3qtrreg] = midptreg(iidx,Fmid);
    
    % Collocation equations
    [Ji,dFdpar_i]=jacobianfile(xreg(1),yreg(:,1),FcnArgs{:});

    for i = 1:Nint
        hi = hreg(i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        [Jip1, dFdpar_ip1]=jacobianfile(xip1,yip1,FcnArgs{:});

        %the interior points
        [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
        [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

        [Jip025, dFdpar_ip025]=jacobianfile(xip025,yip025,FcnArgs{:});
        [Jip05, dFdpar_ip05]=jacobianfile(xip05,yip05,FcnArgs{:});
        [Jip075, dFdpar_ip075]=jacobianfile(xip075,yip075,FcnArgs{:});

        Jip05Jip025=Jip05*Jip025;
        Jip05Jip075=Jip05*Jip075;
        % assembly
        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*numode;
        calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numode;
        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*numode;
        calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);

        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*OCMATCONT.HE.numparameter;
        - (hi/90*( 7*(dFdpar_i+dFdpar_ip1)+32*(dFdpar_ip025+dFdpar_ip075)+12*dFdpar_ip05 )     + ...
            hi*hi/180*( Jip025*(9*dFdpar_i-3*dFdpar_ip1) + Jip075*(3*dFdpar_i-9*dFdpar_ip1) + ...
            Jip05*( 16*(dFdpar_ip025-dFdpar_ip075)+5*(dFdpar_ip1-dFdpar_i) )  )     + ...
            hi*hi*hi/240*( Jip05Jip025*(3*dFdpar_i-dFdpar_ip1)+Jip05Jip075*(3*dFdpar_ip1-dFdpar_i) )  ...
            );
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+odecoord.',1,OCMATCONT.HE.numparameter);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(OCMATCONT.HE.coeffcoordmp+(1:OCMATCONT.HE.numparameter),numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        ROWcounter=ROWcounter+numode;   % next equation

        Ji = Jip1;
        dFdpar_i = dFdpar_ip1;
    end
end
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
J=J(1:gen_c);
% J=sparse(ROW,COL,J,OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables);

%--------------------------------------------------------------------------

function J=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

J   = - In -             ( ...
    hi/90*       ( 7*Ji+27*Jip025+6*Jip05+5*Jip075   ) + ...
    hi*hi/360*   ( 27*Jip05Jip025-5*Jip05Jip075)       + ...
    ( hi*hi/360*(18*Jip025-10*Jip05+6*Jip075) + ...
    hi*hi*hi/240*(3*Jip05Jip025-Jip05Jip075) )*Ji      ...
    );
%--------------------------------------------------------------------------

function J=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

J   = In -               (...
    hi/90*       ( 5*Jip025+6*Jip05+27*Jip075+7*Jip1  ) + ...
    hi*hi/360*   ( 5*Jip05Jip025-27*Jip05Jip075 )       - ...
    ( hi*hi/360*(6*Jip025-10*Jip05+18*Jip075) + ...
    hi*hi*hi/240*(Jip05Jip025-3*Jip05Jip075)  )*Jip1   ...
    );
%--------------------------------------------------------------------------
%---------------------------------------------------------------------------

function [dBCdya,dBCdyb,nCalls,dBCdpar] = BCnumjac(bc,ya,yb,n,npar, ...
    nExtraArgs,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.

% Do not pass info about singular BVPs in ExtraArgs to BC function.
bcArgs = {ya(:),yb(:),n,bc,ExtraArgs{1:nExtraArgs}};
dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
dBCoptions.fac = [];
dBCoptions.vectvars = []; % BC functions not vectorized

bcVal  = bcaux(bcArgs{:});
nCalls = 1;

dBCoptions.diffvar = 1;
[dBCdya,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
nCalls = nCalls + nbc;
dBCoptions.diffvar = 2;
[dBCdyb,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
nCalls = nCalls + nbc;
if npar > 0
    bcArgs = {ya,yb,ExtraArgs{1:nExtraArgs}};
    dBCoptions.thresh = repmat(1e-6,npar,1);
    dBCoptions.diffvar = 3;
    [dBCdpar,ignored,ignored1,nbc] = odenumjac(bc,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;
end

%---------------------------------------------------------------------------

function [dFdy,dFdy_fac,dFdp,dFdp_fac,nFcalls] = Fnumjac(ode,odeArgs,odeVal,...
    Joptions,dPoptions)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.

[dFdy,dFdy_fac,ignored,dFdy_nfcn] = odenumjac(ode,odeArgs,odeVal,Joptions);

[dFdp,dFdp_fac,ignored,dFdp_nfcn] = odenumjac(ode,odeArgs,odeVal,dPoptions);

nFcalls = dFdy_nfcn + dFdp_nfcn;
%--------------------------------------------------------------------------

function res = bcaux(Ya,Yb,n,bcfun,varargin)
ya = reshape(Ya,n,[]);
yb = reshape(Yb,n,[]);
res = feval(bcfun,ya,yb,varargin{:});


%--------------------------------------------------------------------------

function [a025, a05, a075] = midptreg(iidx,mid)
a025 = mid(:,iidx,1);
a05 = mid(:,iidx,2);
a075 = mid(:,iidx,3);

function [a025, a05, a075] = midpti(i, b025, b05, b075)
a025 = b025(:,i);
a05 = b05(:,i);
a075 = b075(:,i);
