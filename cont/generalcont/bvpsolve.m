function sout=bvpsolve(varargin)
%
% BVPCONT main OCMat continuation file for BVPs
%
% BVPCONT(PROBLEMTYPE,INITSOL) argument PROBLEMTYPE is a string
% characterizing the type of problem to be solved:
%   'extremal2ep'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'extremalp2ep'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'extremalt2ep'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'indifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitextremal' ... continuation of a limitpoint solution
%
% INITSOL is an initial function structure to start the continuation,
% returned by an initialization function, e.g. initocmat_AE_EP,
% initocmat_AE_IS (see OCMat manual)
%
% BVPCONT(PROBLEMTYPE,INITSOL,INITTANGENT) INITTANGENT is an initial
% tangent, usually this argument is empty.
%
% BVPCONT(PROBLEMTYPE,INITSOL,INITTANGENT,OPT) the option structure OPT provides a
% multitude of settings for the continuation process and BVP solver (see
% OCMat manual).
%
% SOUT=BVPCONT(...) the structure array consists at least of two elements
% from the intial and last step of the continuation.
%
% The main structure of the bvpcont file is taken from the MatCont file 'cont'

global OCMATCONT
[OCMATCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});
delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% function handles of the actual continuation type
problemhandles=feval(OCMATCONT.problem_func);
OCMATCONT.operatoreq=problemhandles{1};
OCMATCONT.frechetder=problemhandles{2};
OCMATCONT.options=problemhandles{3};
OCMATCONT.defaultprocessor=problemhandles{7};
OCMATCONT.testfunc=problemhandles{8};
OCMATCONT.targetvaluefunc=problemhandles{9};
OCMATCONT.probleminit=problemhandles{10};
OCMATCONT.operatorpfrechet=problemhandles{11};
OCMATCONT.plotcontinuation=problemhandles{21};
OCMATCONT.printcontinuation=problemhandles{30};
OCMATCONT.singmat=problemhandles{13};
OCMATCONT.process=problemhandles{14};
OCMATCONT.locate=problemhandles{15};
OCMATCONT.done=problemhandles{16};
OCMATCONT.adapt=problemhandles{17};
try
    OCMATCONT.dataadaptation=problemhandles{18};
catch
    OCMATCONT.dataadaptation=[];
end
OCMATCONT.stateadapt=problemhandles{19};
OCMATCONT.workspaceadapt=problemhandles{20};
OCMATCONT.formatsolution=problemhandles{22};
OCMATCONT.testadmissibility=problemhandles{23};
OCMATCONT.drearr=problemhandles{24};
OCMATCONT.deval=problemhandles{25};
OCMATCONT.saveintermediate=problemhandles{26};
OCMATCONT.domaindiscretization=problemhandles{28};

ode=problemhandles{4}{1};
bc=problemhandles{4}{2};
if length(problemhandles{4})>=3
    icfun=problemhandles{4}{3};
else
    icfun=[];
end
if getocoptions(opt,'SBVPOC','FJacobian')
    odejac=problemhandles{5}{1};
else
    odejac=[];
end
if getocoptions(opt,'SBVPOC','BCJacobian')
    bcjac=problemhandles{5}{2};
else
    bcjac=[];
end
if getocoptions(opt,'SBVPOC','ICJacobian') && ~isempty(icfun)
    icfunjac=problemhandles{5}{3};
else
    icfunjac=[];
end

%opt=OCMATCONT.options(opt);

OCMATCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCMATCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCMATCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCMATCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCMATCONT.OPTIONS.trm=getocoptions(opt,'NEWTON','TRM');
OCMATCONT.OPTIONS.lambdamin=getocoptions(opt,'NEWTON','LambdaMin');
OCMATCONT.OPTIONS.updatejacfactor=getocoptions(opt,'NEWTON','UpdateJacFactor');
OCMATCONT.OPTIONS.switchtoffnfactor=getocoptions(opt,'NEWTON','SwitchToFFNFactor');
OCMATCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCMATCONT.OPTIONS.display=getocoptions(opt,'NEWTON','Display');
OCMATCONT.OPTIONS.log=getocoptions(opt,'NEWTON','Log');
OCMATCONT.OPTIONS.singularthreshold=getocoptions(opt,'NEWTON','SingularThreshold');
OCMATCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));
OCMATCONT.OPTIONS.admissibletol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
OCMATCONT.OPTIONS.zerotimedifftolerance=getocoptions(opt,'GENERAL','ZeroTimeDiffTolerance');

%OCMATCONT.OPTIONS.residualreductionguard=getocoptions(opt,'SBVPOC','residualReductionGuard');  % to prevent mesh oscillations.
OCMATCONT.OPTIONS.meshadaptation=strcmp(getocoptions(opt,'SBVPOC','MeshAdaptation'),'on');
OCMATCONT.OPTIONS.meshadaptabstol=getocoptions(opt,'SBVPOC','MeshAdaptAbsTol');
OCMATCONT.OPTIONS.meshadaptreltol=getocoptions(opt,'SBVPOC','MeshAdaptRelTol');
OCMATCONT.OPTIONS.xyvectorized=strcmpi(getocoptions(opt,'SBVPOC','Vectorized'),'on');
%OCMATCONT.OPTIONS.meshadaptk=getocoptions(opt,'SBVPOC','MeshAdaptK');
%OCMATCONT.OPTIONS.finemesh=getocoptions(opt,'SBVPOC','MeshAdaptFineMesh');
%OCMATCONT.OPTIONS.meshadaptmaxiter=getocoptions(opt,'SBVPOC','MeshAdaptMaxIter');
OCMATCONT.OPTIONS.maxgridnum=getocoptions(opt,'SBVPOC','NMax');
%OCMATCONT.OPTIONS.maxnewpoints=getocoptions(opt,'SBVPOC','MaxNewPts');

OCMATCONT.OPTIONS.totalrelativedistance=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
OCMATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCMATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCMATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCMATCONT.OPTIONS.increment=getocoptions(opt,'OCCONTARG','Increment');
OCMATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');
OCMATCONT.OPTIONS.saveintermediate=strcmp(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
OCMATCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation

OCMATCONT.OPTIONS.messagedisplay=strcmpi(getocoptions(opt,'INIT','MessageDisplay'),'on');
OCMATCONT.OPTIONS.testmodus=strcmpi(getocoptions(opt,'INIT','TestModus'),'on');

OCMATCONT.bvpmethod=getocoptions(opt,'GENERAL','BVPMethod');

% Handle argument functions and additional arguments
[tmesh0,coeff0]=adapt2solversimple(sol0,ode,odejac,bc,bcjac,icfun,icfunjac);

HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
Singularities=getocoptions(opt,'OCCONTARG','Singularities');
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');
PlotCont=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');
IgnoreSings=getocoptions(opt,'OCCONTARG','IgnoreSingularity');
ExitOnTargetValue=strcmpi(getocoptions(opt,'OCCONTARG','ExitOnTargetValue'),'on');
Locators=getocoptions(opt,'OCCONTARG','Locators');
Userfunctions=getocoptions(opt,'OCCONTARG','Userfunctions');
MaxContStep=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
Backward=getocoptions(opt,'OCCONTARG','Backward');
CheckAdmissibility=strcmpi(getocoptions(opt,'OCCONTARG','CheckAdmissibility'),'on');
Adapt=getocoptions(opt,'OCCONTARG','Adapt');
MeshAdaptation=getocoptions(opt,'SBVPOC','MeshAdaptation');

% used options in the continuation file
stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth'); % initialize step width
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

OCMATCONT.nActTest=0;
% initialize user workspace
% if WorkSpace
%     if OCMATCONT.probleminit(tmesh0,coeff0,tangent0)~=0
%         ocmaterror('Initializer failed.');
%     end
% end
newmesh=true;
while newmesh
    [tmesh,coeff]=newtcorr4bvpsolve(tmesh0,coeff0);
    if ~isempty(tmesh)
        [tmesh0,coeff0,newmesh]=AdaptMesh(tmesh,coeff);
        DataAdaptation(tmesh0);
        OCMATCONT.HE.contparametercoord=[];
    else
        break
    end
end

if ~isempty(coeff)
    s.index=1;
    s.label='99';
    s.data.sol=OCMATCONT.formatsolution(tmesh0,coeff0,[]);
    s.msg='This is the last solution of the BVP continuation';
    sout=s;
else
    sout=[];
end

if SaveIntermediate
    bvpout.tmesh=tmesh0;
    bvpout.coeff=coeff0;
    bvpout.tangent=[];
    failed=OCMATCONT.saveintermediate(sout,bvpout,1);
end
end
%----------------------------
%
% Check admissibility of solution
%
%----------------------------

function [failed,infoS,labelS]=EvalAdmissibilityFunc(tmesh,coeff,tangent)
global OCMATCONT

[failed,infoS,labelS]=OCMATCONT.testadmissibility(tmesh,coeff,tangent);
end
%--< END OF EVALADMISSIBILITYFUNC >--


%---------------------------------
%
% Plot the solution
%
%
function b=PlotContinuationProcess(tmesh,coeff,tangent)
global OCMATCONT
b=OCMATCONT.plotcontinuation(tmesh,coeff,tangent);

%-----------------------------------
%
% Command line parser
%
%-----------------------------------
end
function [problemfunc,solinit,opt]=ParseCommandLine(problemfunc,solinit,varargin)

if nargin < 2
    ocmaterror('wrong number of input arguments');
end
opt=[];
if nargin > 2
    solinit.solverinfo.tangent=varargin{1};
    if nargin > 3
        opt=varargin{2};
    end
end
if isempty(opt)
    opt=defaultocoptions;
end
%--< END OF CMDL PARSER >--
end

function [tmesh,coeff]=adapt2solversimple(solinit,ode,odejac,bc,bcjac,icfun,icfunjac)

clear global OCBVP
global OCMATCONT OCBVP

bvpmethod=OCMATCONT.bvpmethod;
% odefinal=ode;
% bcfinal=bc;
% icfinal=ic;
% jacfinal=odejac;
% bcjacfinal=bcjac;
% icjacfinal=icjac;

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

OCMATCONT.HE.numarc=numel(solinit.arcarg);
OCMATCONT.HE.arcindex=arcarg2arcindex(solinit.arcarg);
OCMATCONT.HE.arcarg=solinit.arcarg;

% edge [arcarg1 arcarg2] denotes the switch from arc with
% arc argument arcarg1 to arc argument arcarg2, depending on this
% information e.g. the boundary condtions have to be chosen
OCMATCONT.HE.edge=[OCMATCONT.HE.arcarg(1:end-1);OCMATCONT.HE.arcarg(2:end)];
OCMATCONT.HE.numparameter=numel(solinit.parameters);%additional/continuation parameters
OCBVP.npar=OCMATCONT.HE.numparameter; % nmuber of free paramteres
OCMATCONT.codimension=0;

% set time mesh
switch bvpmethod
    case {'bvp4c','bvp6c','gbvp4c'}
        numode=size(solinit.y,1);
        tmesh=solinit.x;
        coeff=solinit.y(:);
        coeff=[coeff;solinit.parameters(:)];
        OCMATCONT.HE.numdvariables=numel(coeff);
        OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;
end

OCMATCONT.HE.TIMEDDATA.nummesh=numel(tmesh);
OCMATCONT.HE.TIMEDDATA.diffmesh=diff(tmesh);
OCMATCONT.HE.TIMEDDATA.mesh=tmesh;

arcposition=find(diff(solinit.x)==0);
OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 arcposition+1];
OCMATCONT.HE.TIMEDDATA.rightarcindex=[arcposition OCMATCONT.HE.TIMEDDATA.nummesh];
OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;

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
        OCBVP.nBCs=OCMATCONT.HE.totalnumboundarycondition;%-1;%-OCBVP.numic;

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
    jacfinal(t(1),y(:,1),1,freepar,modelpar);
catch
    jacfinal=[];
    OCBVP.multiarccalc=multiarccalc;
end
try
    multiarccalc=OCBVP.multiarccalc;
    bcjacfinal(y(:,OCBVP.Lidx),y(:,OCBVP.Ridx),freepar,modelpar);
catch
    bcjacfinal=[];
    OCBVP.multiarccalc=multiarccalc;
end

% options for numerical differentiation
Joptions=[];
dPoptions=[];
dBCoptions=[];
if 1%isempty(jacfinal)
    Joptions.diffvar=2;  % dF(x,y)/dy
    Joptions.fac = [];
    if OCMATCONT.OPTIONS.xyvectorized
        Joptions.vectvars=[1,2];
    else
        Joptions.vectvars=[];
    end
    dPoptions.diffvar=4;  % dF(x,y,region,p)/dp
    dPoptions.vectvars=[]; % no vectorization for parameter Jacobian
    dPoptions.fac = [];
end
if 1%isempty(bcjacfinal)
    dBCoptions.vectvars=[];
    dBCoptions.fac_dya=[];
    dBCoptions.fac_dyb=[];
end
OCBVP.Joptions=Joptions;
OCBVP.dPoptions=dPoptions;
OCBVP.dBCoptions=dBCoptions;
OCBVP.averageJac=isempty(jacfinal);
OCBVP.nparmc=0;

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

OCMATCONT.ode=odefinal;
OCMATCONT.bc=bcfinal;
OCMATCONT.icfun=icfunfinal;
OCMATCONT.odejac=jacfinal;
OCMATCONT.bcjac=bcjacfinal;
OCMATCONT.icfunjac=icfunjacfinal;

% ---------------------------------------------------------
% Nested functions
% ---------------------------------------------------------

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
        OCMATCONT.HE.contparametercoord=[];

        if OCMATCONT.HE.numarc>1 && ~isempty(arcposition)
            OCBVP.multipointbvp=true;
            OCBVP.numarc=OCMATCONT.HE.numarc;
            OCBVP.multiarccalc=true;
        else
            OCBVP.multipointbvp=false;
            OCBVP.numarc=1;
            OCBVP.multiarccalc=true;
        end
        OCMATCONT.HE.totalnumboundarycondition=OCMATCONT.HE.numarc*numode+OCMATCONT.HE.numparameter;
        OCBVP.nBCs=OCMATCONT.HE.totalnumboundarycondition;
        OCBVP.Lidx=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        OCBVP.Ridx=OCMATCONT.HE.TIMEDDATA.rightarcindex;
        OCBVP.Nint=OCMATCONT.HE.TIMEDDATA.nummeshintv;
        OCBVP.rows=OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.numode;   % define the action area
        OCBVP.arcposition=arcposition;
    end
    function f=odeIntegralConstraint(x,y,arcarg,p,par)
        % add integral constraint as a further ODE
        f=[ode(x,y(1:numodeexcnumic,:),arcarg,p,par);
            sum(icfun(x,arcarg,p,par).*y(1:numodeexcnumic,:))];
    end  % odeIntegralConstraint

% ---------------------------------------------------------

    function res=bcIntegralConstraint(ya,yb,p,par)
        res=[bc(ya(1:numodeexcnumic),yb(1:numodeexcnumic),p,par);
            ya(numodeexcnumic+1:numode);
            yb(numodeexcnumic+1:numode)];
    end  % bcIntegralConstraint

% ---------------------------------------------------------

    function [J Jpar]=jacIntegralConstraint(x,y,arcarg,p,par)
        [dfdy,dfdp]=odejac(x,y(1:numodeexcnumic,:),arcarg,p,par);
        [dicdy,dicdp]=icfunjac(x,y,arcarg,p,par);
        % add Jacobian for the integral constraint
        J=[dfdy zeros(numodeexcnumic,numic);dicdy zeros(numic)];
        Jpar=[dfdp;dicdp];
    end  % jacIntegralConstraint

% ---------------------------------------------------------

    function [dya,dyb,dypar]=bcjacIntegralConstraint(ya,yb,p,par)
        [dbcdya,dbcdyb,dbcdp]=bcjac(ya(1:numodeexcnumic,:),yb(1:numodeexcnumic,:),p,par);
        dya=[dbcdya;ones(numic,1)];
        dyb=[dbcdyb;ones(numic,1)];
        dypar=[dbcdp;zeros(numic,1)];
    end  % bcjacIntegralConstraint

% ---------------------------------------------------------

    function f=odeParameters(x,y,arcarg,p,par)
        p=y(numodeexcnumic+1:numodeexcnumic+nparam,1);  % extract p from y(:,1)
        % add trivial equations for unknown parameters
        f=[ ode(x,y(1:numodeexcnumic,:),arcarg,p,par);
            zeros(nparam,numel(x))];
    end  % odeParameters

% ---------------------------------------------------------

    function res=bcParameters(ya,yb,p,par)
        p=ya(numode+1:numode+nparam);   % extract p from ya
        res=bc(ya(1:numode),yb(1:numode),p,par);
    end  % bcParameters

% ---------------------------------------------------------

    function J=jacParameters(x,y,arcarg,p,par)
        [dfdy,dfdp]=odejac(x,y(1:numode,:),arcarg,p,par);
        % add trivial equations for unknown parameters
        J=[dfdy, dfdp;
            zeros(nparam,numode+nparam)];
    end  % jacParameters

% ---------------------------------------------------------

    function [dya,dyb]=bcjacParameters(ya,yb,p,par)
        [dbcdya,dbcdyb,dbcdp]=bcjac(ya(1:numode,:),yb(1:numode,:),p,par);
        dya=[dbcdya, dbcdp];
        dyb=[dbcdyb, zeros(numode,nparam)];
    end  % bcjacParameters

% ---------------------------------------------------------

    function res=bcCollocation(tmesh,ya,yb,za,zb,p,par)
        [yl,yr]=boundarycollocationpoints(tmesh,ya,yb,za,zb);   % extract p from ya
        res=bc(yl,yr,p,par);
    end  % bcParameters

% ---------------------------------------------------------

    function [dya,dyb,dypar]=bcjacCollocation(tmesh,ya,yb,za,zb,p,par)
        [yl,yr]=boundarycollocationpoints(tmesh,ya,yb,za,zb);   % extract p from ya
        [dya,dyb,dypar]=bcjac(yl,yr,p,par);
    end  % bcParameters

end




%------------------------------------------------
%
% Call method specific mesh adaptation
%
%------------------------------------------------


function [tmesh,coeff,newmesh]=AdaptMesh(tmesh,coeff)
global OCMATCONT
newmesh=false;
[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            return
        end
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if verifyresidual_bvp5c(max(res))
            return
        end
    case 'bvp4c'
        res=residual_bvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            return
        end
    case 'gbvp4c'
        res=residual_gbvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            return
        end
end
% Detect mesh oscillations:  Was there a mesh with
% the same number of nodes and a similar residual?
% residualReduction = abs(OCMATCONT.meshHistory((OCMATCONT.meshHistory(:,1) == OCMATCONT.HE.TIMEDDATA.nummesh),2) - maxres)/maxres;
oscLikely = false;%any( residualReduction < OCMATCONT.OPTIONS.residualreductionguard);

switch OCMATCONT.bvpmethod
    case 'bvp6c'
        [tmesh,ynew]=meshadaptation_bvp6c(t,y,[],res,oscLikely);
        coeff=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tmesh,ynew]=meshadaptation_bvp5c(t,y,[],res,oscLikely,freepar,modelpar,OCMATCONT.ode);
        coeff=[ynew(:)];
    case 'bvp4c'
        [tmesh,ynew]=meshadaptation_bvp4c(t,y,[],res,oscLikely);
        coeff=[ynew(:);freepar(:)];
    case 'gbvp4c'
        [tmesh,ynew,tangent]=meshadaptation_gbvp4c(t,y,[],res,~oscLikely);
        dataadaptation_gbvp4c(tmesh)
        coeff=ynew(OCMATCONT.HE.DDATA.meshvalcoord);
        coeff=[coeff(:);freepar(:)];
end

newmesh=true;
end