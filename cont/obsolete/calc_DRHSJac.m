function Jac=calc_DRHSJac(x,y,freepar,modelpar,map,bc,ic,mapjac,bcjac,icjac)

global OCBVP

FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.

threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(OCBVP.nummap,1));
OCBVP.dPoptions.thresh = threshval(ones(OCBVP.npar,1));

F=OCBVP.F;

% BC points
ya = y(:,OCBVP.Lidx);
yb = y(:,OCBVP.Ridx);

rows = OCBVP.rows;   % define the action area
cols = OCBVP.cols;   % in the global Jacobian
OCBVP.nNm1=OCBVP.nN;
Jac = spalloc(OCBVP.nNm1+OCBVP.nparmc,OCBVP.nNm1+OCBVP.npar,(OCBVP.nNm1+OCBVP.npar)*(2*OCBVP.nummap+OCBVP.npar));  % sparse storage
last_cols = zeros(OCBVP.nNm1+OCBVP.nparmc,OCBVP.npar);   % accumulator

if isempty(bcjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = BCnumjac(bc,ya,yb,FcnArgs(2:3));
else  % use analytical Jacobian
    [dGdya,dGdyb,dGdpar] = bcjac(ya,yb,FcnArgs{2:3});
end
last_cols(1:OCBVP.nBCs,:) = dGdpar;

if OCBVP.sumconstraint
    if isempty(icjac)  % use numerical approx
        [JICdy,JICdp]=ICnumjac(ic,{x,y,FcnArgs{2:end}});
    else
        [JICdy,JICdp]=icjac(x,y,FcnArgs{2:end});
    end
    Jac(OCBVP.nNm1+OCBVP.nparmc-OCBVP.nICs+1:OCBVP.nNm1+OCBVP.nparmc,1:OCBVP.nNm1+OCBVP.nparmc)=JICdy;
    last_cols(OCBVP.nNm1+OCBVP.nparmc-OCBVP.nICs+1:OCBVP.nNm1+OCBVP.nparmc,:)=JICdp;
end

% Collocation equations
for arc = 1:OCBVP.numarc

    % Left BC
    Jac(1:OCBVP.nBCs,cols) = dGdya(:,(arc-1)*OCBVP.nummap+(1:OCBVP.nummap));

    FcnArgs{1} = arc;

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    xreg = x(xidx(1:end-1));
    yreg = y(:,xidx);
    Freg = F(:,xidx(1:end-1));

    % Collocation equations
    if isempty(mapjac)  % use numerical approx
        for ii = 1:OCBVP.Nint(arc)
            % the left mesh point
            OCBVP.Fref=Freg(:,ii);
            OCBVP.Joptions.diffvar=2;  % dF(x,y)/dy
            OCBVP.Joptions.vectvars=[];
            Ji=Fnumjac(map,{xreg(ii+1),yreg(:,[ii ii+1]),FcnArgs{:}});
            OCBVP.Joptions.diffvar=3;  % dF(x,y,p)/dp
            OCBVP.Joptions.vectvars=[];
            [Jip1,dFdpar_ip1]=Fnumjac(map,{ xreg(ii+1),yreg(:,[ii ii+1]),FcnArgs{:}});

            % assembly
            Jac(rows,cols) =Ji;
            cols = cols + OCBVP.nummap;
            Jac(rows,cols) =Jip1;
            last_cols(rows,:) =dFdpar_ip1;
            rows = rows+OCBVP.nummap;   % next equation
        end

    else % use analytical Jacobian
        for ii = 1:OCBVP.Nint(arc)
            % the left mesh point
            [Ji dFdpar_ip1]=mapjac(xreg(ii+1),yreg(:,[ii ii+1]),FcnArgs{:});

            Jac(rows,[cols cols + OCBVP.nummap]) =Ji;
            cols = cols + OCBVP.nummap;
            last_cols(rows,:) =dFdpar_ip1;
            rows = rows+OCBVP.nummap;   % next equation
        end
    end
    % Right BC
    Jac(1:OCBVP.nBCs,cols) = dGdyb(:,(arc-1)*OCBVP.nummap+(1:OCBVP.nummap));
    %cols = cols + OCBVP.nummap;
end
Jac(:,end-OCBVP.npar+1:end) = last_cols;  % accumulated

function [dBCdya,dBCdyb,nCalls,dBCdpar]=BCnumjac(bc,ya,yb,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.
global OCBVP
dBCoptions=OCBVP.dBCoptions;

if OCBVP.multipointbvp
    % make ya and yb columns vectors
    bcArgs={ya(:),yb(:),OCBVP.neqn,bc,ExtraArgs{:}};
    dBCoptions.diffvar=1;
    [dBCdya,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nbc;
    dBCoptions.diffvar=2;
    [dBCdyb,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
else
    bcArgs={ya,yb,ExtraArgs{:}};
    dBCoptions.diffvar=1;
    [dBCdya,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nbc;
    dBCoptions.diffvar=2;
    [dBCdyb,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
end
if OCBVP.explicitparameterdependence
    bcArgs={ya,yb,ExtraArgs{:}};
    dBCoptions.diffvar=3;
    [dBCdpar,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
else
    dBCdpar=[];
end

function [dFdy,dFdp]=Fnumjac(map,ExtraArgs)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
global OCBVP
mapArgs={ExtraArgs{1},ExtraArgs{2}(:,1),ExtraArgs{2}(:,2),map,ExtraArgs{3:end}};
dFdy=numjaccsd(@mapaux,mapArgs,OCBVP.neqn,OCBVP.Joptions);
if OCBVP.explicitparameterdependence
   dFdp=numjaccsd(@mapaux,mapArgs,OCBVP.neqn,OCBVP.dPoptions);
else
    dFdp=[];
end

function [dICdy,dICdp]=ICnumjac(ic,ExtraArgs)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
global OCBVP
dICdp=[];
mapArgs={ExtraArgs{1},ExtraArgs{2}(:),OCBVP.neqn,ic,ExtraArgs{3:end}};
dICdy=numjaccsd(@icaux,mapArgs,OCBVP.nICs,OCBVP.dICoptions);
if OCBVP.explicitparameterdependence
    mapArgs={ExtraArgs{1},ExtraArgs{2}(:),OCBVP.neqn,ic,ExtraArgs{3:end}};
    dICdp=numjaccsd(@icaux,mapArgs,OCBVP.nICs,OCBVP.dICdPoptions);  
end

%---------------------------------------------------------------------------
function res=mapaux(x,yi,yip1,mapfun,varargin)
res=mapfun(x,[yi yip1],varargin{:});

%---------------------------------------------------------------------------
function res=icaux(x,y,n,icfun,varargin)
y=reshape(y,n,[]);
res=icfun(x,y,varargin{:});

%---------------------------------------------------------------------------
function res=bcaux(Ya,Yb,n,bcfun,varargin)
ya=reshape(Ya,n,[]);
yb=reshape(Yb,n,[]);
res=bcfun(ya,yb,varargin{:});
