function [tmesh,coeff,tangent]=adapt2solver(solinit,ode,odejac,bc,bcjac,icfun,icfunjac,odehess,bchess,odetensor3,bctensor3)

clear global OCBVP
global OCMATCONT OCBVP

bvpmethod=OCMATCONT.bvpmethod;
% odefinal=ode;
% bcfinal=bc;
% icfinal=ic;
% jacfinal=odejac;
% bcjacfinal=bcjac;
% icjacfinal=icjac;
odehessfinal=[];
bchessfinal=[];
odetensor3final=[];
bctensor3final=[];
% Validate arguments
if ~isfield(solinit,'x')
    msg=sprintf('The field ''x'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
elseif ~isfield(solinit,'y')
    msg=sprintf('The field ''y'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
elseif ~isfield(solinit,'parameters')
    msg=sprintf('The field ''parameters'' not present in SOLINIT.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

if length(solinit.x) < 2
    msg=sprintf('SOLINIT.x must contain at least the two end points.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

if any(sign(solinit.x(end)-solinit.x(1)) * diff(solinit.x) < 0)
    msg=sprintf('The entries in SOLINIT.x must increase or decrease.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

if isempty(solinit.y)
    msg=sprintf('No initial guess provided in SOLINIT.y.');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end
if size(solinit.y,2) ~= length(solinit.x)
    msg=sprintf('SOLINIT.y not consistent with SOLINIT.x.');
    error('MATLAB:bvparguments:SolXSolYSizeMismatch','%s  %s',msg);
end

if isfield(solinit,'solverinfo') && isfield(solinit.solverinfo,'tangent')
    tangent=solinit.solverinfo.tangent;
else
    tangent=[];
end
if ~isfield(OCMATCONT,'TargetValueNum')
    OCMATCONT.TargetValueNum=1;
end
OCMATCONT.HE.numarc=numel(solinit.arcarg);
OCMATCONT.HE.arcindex=arcarg2arcindex(solinit.arcarg);
OCMATCONT.HE.arcarg=solinit.arcarg;

% edge [arcarg1 arcarg2] denotes the switch from arc with
% arc argument arcarg1 to arc argument arcarg2, depending on this
% information e.g. the boundary condtions have to be chosen
OCMATCONT.HE.edge=[OCMATCONT.HE.arcarg(1:end-1);OCMATCONT.HE.arcarg(2:end)];
OCMATCONT.HE.numparameter=numel(solinit.parameters);%additional/continuation parameters
OCMATCONT.HE.numparametermc=OCMATCONT.HE.numparameter-1;%OCMATCONT.codimension;%additional parameters minus continuation parameter
OCBVP.npar=OCMATCONT.HE.numparameter; % nmuber of free paramteres
OCBVP.nparmc=OCMATCONT.HE.numparameter-1;%-OCMATCONT.codimension; % nmuber of free paramteres excluded the number of continuation parameters
OCBVP.nparmcod=OCMATCONT.HE.numparameter-OCMATCONT.codimension; % nmuber of free paramteres excluded the number of continuation parameters

% set time mesh
numode=size(solinit.y,1);
tmesh=solinit.x;
coeff=solinit.y(:);
coeff=[coeff;solinit.parameters(:)];
OCMATCONT.HE.numdvariables=numel(coeff);
OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
if ~isempty(tangent) && (length(coeff) ~= length(tangent))
    msg=sprintf('tangent vector and coefficient are not consistent');
    ocmaterror('MATLAB:bvpcont:NoXInSolinit','%s',msg);
end

OCMATCONT.HE.TIMEDDATA.nummesh=numel(tmesh);
OCMATCONT.HE.TIMEDDATA.diffmesh=diff(tmesh);
OCMATCONT.HE.TIMEDDATA.mesh=tmesh;

arcposition=find(diff(solinit.x)==0);
OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 arcposition+1];
OCMATCONT.HE.TIMEDDATA.rightarcindex=[arcposition OCMATCONT.HE.TIMEDDATA.nummesh];
OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;


% Function handles
odefinal=@(x,y,region,p,par)ode(x,y,region,p,par);
if isa(odejac,'function_handle')
    jacfinal=@(x,y,region,p,par)odejac(x,y,region,p,par);
end
bcfinal=@(ya,yb,p,par) bc(ya,yb,p,par);
if isa(bcjac,'function_handle')
    bcjacfinal=@(ya,yb,p,par) bcjac(ya,yb,p,par);
end
if ~isempty(icfun)
    numic=OCMATCONT.numintegralconstraint;
    OCBVP.numic=numic;
    % handle integral constraints of the form
    % int_0^T(icfun(x,region,par)*y(x))=0
    icfunfinal=@(x,y,region,p,par)icfun(x,y,region,p,par);
    if isa(icfunjac,'function_handle')
        icfunjacfinal=@(x,y,region,p,par)icfunjac(x,y,region,p,par);
    else
        icfunjacfinal=[];
    end
else
    OCBVP.numic=0;
    icfunfinal=[];
    icfunjacfinal=[];
end
if isa(odehess,'function_handle')
    odehessfinal=@(x,y,region,p,par)odehess(x,y,region,p,par);
end
if isa(bchess,'function_handle')
    bchessfinal=@(ya,yb,p,par)bchess(ya,yb,p,par);
end

switch bvpmethod
    case {'bvp4c','bvp6c'}
        initbvp4c()
    case 'gbvp4c'
        initbvp4c()
        OCBVP.numode=OCMATCONT.numode;
        OCBVP.maxnumode=OCMATCONT.maxnumode;
        
        OCBVP.Nint=zeros(1,size(OCBVP.Lidx,2));
        OCBVP.In=cell(1,size(OCBVP.Lidx,2));
        counter=0;
        counter1=0;
        totidx=[];
        keepbccolumns=[];
        for ii=1:size(OCBVP.Lidx,2)
            Y=zeros(OCBVP.maxnumode,OCBVP.Ridx(ii)-OCBVP.Lidx(ii)+1);
            Y(1:OCBVP.numode(ii),:)=1;
            idx=find(Y);
            totidx=[totidx;idx+counter];
            keepbccolumns=[keepbccolumns counter1+(1:OCBVP.numode(ii))];
            OCBVP.Nint(ii)=OCBVP.Ridx(ii)-OCBVP.Lidx(ii);
            OCBVP.In{ii}=eye(OCBVP.numode(ii));
            counter=counter+OCBVP.maxnumode*(OCBVP.Nint(ii)+1);
            counter1=counter1+OCBVP.maxnumode;
        end
        OCBVP.keepbccolumns=keepbccolumns;
        OCMATCONT.HE.DDATA.meshvalcoord=totidx;
        coeff=solinit.y(totidx);
        OCMATCONT.HE.parametercoord=length(coeff);
        OCMATCONT.HE.ycoord=1:OCMATCONT.HE.parametercoord;
        coeff=[coeff;solinit.parameters(:)];
        OCMATCONT.HE.parametercoord=OCMATCONT.HE.parametercoord+(1:length(solinit.parameters));
        OCBVP.nN=sum((OCBVP.Nint+1).*OCBVP.numode);
        OCMATCONT.HE.totalnumboundarycondition=sum(OCBVP.numode)+OCMATCONT.HE.numparameter-OCMATCONT.codimension;
        OCBVP.nBCs=OCMATCONT.HE.totalnumboundarycondition-OCBVP.numic;
        
        OCBVP.cols=1:OCBVP.numode(1);             % in the global Jacobian
        OCBVP.rows=OCBVP.nBCs+(1:OCBVP.numode(1));
        OCMATCONT.HE.numdvariables=sum((OCBVP.Nint+1).*OCBVP.numode)+OCMATCONT.HE.numparameter;
        OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
        OCMATCONT.HE.contparametercoord=length(coeff);

end
if ~OCBVP.multipointbvp
    if ~isempty(arcposition)
        ocmaterror('For two-point BVPs the timemesh has to be strictly de/increasing.')
    end
    rowcounter=0;
    rowarcposition=zeros(2,OCMATCONT.HE.numarc);
    for ii=1:OCMATCONT.HE.numarc
        rowcounter_start=rowcounter+1;
        rowcounter=rowcounter+numode;
        rowarcposition(1:2,ii)=[rowcounter_start;rowcounter];
    end
    OCMATCONT.HE.arcrowindex=rowarcposition;
    OCMATCONT.HE.bcindex=reshape(1:numode*OCMATCONT.HE.numarc,numode,OCMATCONT.HE.numarc);
end
warnoffId={'MATLAB:singularMatrix','MATLAB:nearlySingularMatrix'};
for ii=1:length(warnoffId)
    warnstat(ii)=warning('query',warnoffId{ii});
    warnoff(ii)=warnstat(ii);
    warnoff(ii).state='off';
end
OCBVP.warnstat=warnstat;
OCBVP.warnoff=warnoff;
OCBVP.warnoffId=warnoffId;

% test functions
[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
try
    multiarccalc=OCBVP.multiarccalc;
    jacfinal(t(1),y(1:OCBVP.numode(1),1),1,freepar,modelpar);
catch
    jacfinal=[];
    OCBVP.multiarccalc=multiarccalc;
end
try
    multiarccalc=OCBVP.multiarccalc;
    bcjacfinal(y(1:OCBVP.numode(1),OCBVP.Lidx),y(1:OCBVP.numode(1),OCBVP.Ridx),freepar,modelpar);
catch
    bcjacfinal=[];
    OCBVP.multiarccalc=multiarccalc;
end

% options for numerical differentiation
Joptions=[];
dPoptions=[];
dBCoptions=[];
Joptions.diffvar=2;  % dF(x,y)/dy
Joptions.fac = [];
Joptions.thresh=1e-6;
if OCMATCONT.OPTIONS.xyvectorized
    Joptions.vectvars=[1,2];
else
    Joptions.vectvars=[];
end
dPoptions.diffvar=4;  % dF(x,y,region,p)/dp
dPoptions.vectvars=[]; % no vectorization for parameter Jacobian
dPoptions.fac = [];
dBCoptions.vectvars=[];
dBCoptions.fac_dya=[];
dBCoptions.fac_dyb=[];
OCBVP.Joptions=Joptions;
OCBVP.dPoptions=dPoptions;
OCBVP.dBCoptions=dBCoptions;
OCBVP.averageJac=isempty(jacfinal);

OCMATCONT.ode=odefinal;
OCMATCONT.bc=bcfinal;
OCMATCONT.icfun=icfunfinal;
OCMATCONT.odejac=jacfinal;
OCMATCONT.bcjac=bcjacfinal;
OCMATCONT.icfunjac=icfunjacfinal;
OCMATCONT.odehess=odehessfinal;
OCMATCONT.bchess=bchessfinal;
OCMATCONT.odetensor3=odetensor3final;
OCMATCONT.bctensor3=bctensor3final;


    function initbvp4c()


        %OCBVP.MaxNewPts=OCMATCONT.OPTIONS.maxnewpoints;
        OCBVP.numode=numode;
        OCBVP.N=numel(tmesh);
        OCBVP.neqn=OCBVP.numode;
        OCBVP.explicitparameterdependence=true;
        OCBVP.nN=OCBVP.numode*OCBVP.N;
        OCBVP.In=eye(OCBVP.numode);
        OCBVP.cols=1:OCBVP.numode;             % in the global Jacobian
        OCBVP.rtol=OCMATCONT.OPTIONS.meshadaptreltol;
        OCBVP.Nmax=OCMATCONT.OPTIONS.maxgridnum;
        OCBVP.threshold=OCMATCONT.OPTIONS.meshadaptabstol/OCMATCONT.OPTIONS.meshadaptreltol;
        if isscalar(OCBVP.threshold)
            OCBVP.threshold=OCBVP.threshold(ones(OCBVP.neqn,1));
        end

        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCBVP.N*OCBVP.neqn,OCBVP.neqn,OCBVP.N);
        OCMATCONT.HE.parametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:OCMATCONT.HE.numparameter).';
        OCMATCONT.HE.parametermccoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:OCMATCONT.HE.numparametermc).';
        OCMATCONT.HE.parametermcodcoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:OCMATCONT.HE.numparameter-OCMATCONT.codimension).';
        OCMATCONT.HE.contparametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+OCMATCONT.HE.numparameter;%-1+(1:OCMATCONT.codimension);

        if OCMATCONT.HE.numarc>1 && ~isempty(arcposition)
            OCBVP.multipointbvp=true;
            OCBVP.numarc=OCMATCONT.HE.numarc;
            OCBVP.multiarccalc=true;
        else
            OCBVP.multipointbvp=false;
            OCBVP.numarc=1;
            OCBVP.multiarccalc=true;
        end
        OCMATCONT.HE.totalnumboundarycondition=OCMATCONT.HE.numarc*numode+OCMATCONT.HE.numparameter-OCMATCONT.codimension;
        OCBVP.nBCs=OCMATCONT.HE.totalnumboundarycondition-OCBVP.numic;
        OCBVP.Lidx=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        OCBVP.Ridx=OCMATCONT.HE.TIMEDDATA.rightarcindex;
        OCBVP.Nint=OCMATCONT.HE.TIMEDDATA.nummeshintv;
        OCBVP.rows=OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numode;   % define the action area
        OCBVP.arcposition=arcposition;
    end
end
