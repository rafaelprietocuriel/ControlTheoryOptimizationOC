function sol=initocmat_AE_EP(ocObj,ocEP,contcoordinate,targetvalue,opt,varargin) 
%
% INITOCMAT_AE_EP initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_EP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
% saddle path calculation is initialized. 
% OCOBJ          ... corresponding optimal control model
% OCEP           ... equilibrium point (dynprimitive) (hat-x) with (local)
%                    stable manifold of dimension k 
% CONTCOORDINATE ... coordinates i_1,...,i_k of the continuation variable
%                    (usually state coordinate(s) in optimal control
%                    problems)  
% TARGETVALUE    ... determines direction of the continuation (x_j^0,
%                    j=i_1,...,i_k) 
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields 
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and OCMATAE are
% initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SHORT EXPLANATION OF THE MATHEMATICAL BACKGROUND
%
% The underlying mathematical problem formulation is to find a trajectory
% (time path) x(t)=(x_1(t),...,x_N(t)) with initial condition x_j(0)=x_j^0,
% j=i_1,...,i_k and lim_{t\to\infty}x(t)=hat-x (convergence to the
% equilibrium). 
%
% The problem is solved using a continuation algorithm, where the
% continuation is done for the initial condition
%   x_j(0)=x_j^0, j=i_1,...,i_k
% and the continuation parameter 'mu' is defined as
%   x_j(0)=x_j^0*mu+(1-mu)*hat-x_j, j=i_1,...,i_k.
% Thus, for mu=0 we have
%   x_j(0)=hat-x_j, j=i_1,...,i_k
% and for mu=1
%   x_j(0)=x_j^0
% The end condition, convergence to the equilibrium, is reformulated in a
% way that allows a numerical treatment. The default way is the truncation
% of the infinite time to a finite time 'T' and the condtion that the end
% point x(T) ends on the linearized stable manifold (stable eigenspace).
%
% This means that at the start of the continuation the equilibrium path
% (constant solution at the equilibrium) trivially satisfies the initial
% and end condition. Therefore, with the provision of the equilibrium the
% initial solution is given as well. In that sense the OCEP argument
% performs two tasks. The searched for solution converges to OCEP and OCEP
% is the initial solution (mu=0) of the continuation process.
%
% The denomination as TARGETVALUE maybe misleading but from the
% continuation point of view it denotes the target. From the problem
% perspective it denotes the initial state values of the searched for
% solution. 
% 
% SOL=INITOCMAT_AE_EP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE,OPT) with the
% option structure OPT the threshold 'ZeroDeviationTolerance' and initial
% number of discretization points 'TrivialArcMeshNum' for the equilibrium
% solution can be changed. 
%   ZeroDeviationTolerance ... provides the tolerance to classify an 
%                              eigenvalue numerically as zero.
%   TrivialArcMeshNum      ... provides the number of points for the
%                              constant solution at the equilbrium. 
%
% SOL=INITOCMAT_AE_EP(...,'TruncationTime',T) the truncation of the
% infinite time horizon to the finite time T 
% 
% SOL=INITOCMAT_AE_EP(...,'PathType',p)
%   p='s' (default) stable saddle-path calculation
%   p='u' unstable saddle-path calculation

%clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
if isdifferentialgame(ocObj)
    sol=initocmatdg_AE_EP(ocObj,ocEP,contcoordinate,targetvalue,opt,varargin{:});
    return
end
sol=[];
pathtype='';
asymptoticapproximation=[];
infinity=[];
freevector=[];
freeendtime=[];
fixedcoord=[];
excludecoord=[];
stopcriterion=[];
objectivevaluecalc=[];
stateconstraintflag=[];
asymptoticmatrix=[];
excludecoordinate4ep=[];
exogenousfunction=[];
userbc=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocEP)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if 0%~testconsistency(ocEP,ocObj)
    ocmaterror('The equilibrium and oc model are not consistent.')
end
if nargin==4
    opt=defaultocoptions;
end
asymptoticapproximationidx=find(strcmpi(varargin,'truncationtime'));
infinityidx=find(strcmpi(varargin,'inftimetransformation'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
freevectoridx=find(strcmpi(varargin,'freevector'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
excludecoordidx=find(strcmpi(varargin,'excludecoord'),1);
fixedcoordidx=find(strcmpi(varargin,'fixedcoord'),1);
stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
stateconstraintidx=find(strcmpi(varargin,'stateconstraint'));
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'));
excludecoordinate4epidx=find(strcmpi(varargin,'excludecoordinate4ep'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
userbcidx=find(strcmpi(varargin,'userbc'));

if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(asymptoticapproximationidx)
    asymptoticapproximation=varargin{asymptoticapproximationidx+1};
end
if ~isempty(infinityidx)
    infinity=varargin{infinityidx+1};
end
if ~isempty(infinity) && infinity
    asymptoticapproximation=1;
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(pathtypeidx)
    pathtype=varargin{pathtypeidx+1};
end
if ~isempty(fixedcoordidx)
    fixedcoord=varargin{fixedcoordidx+1};
end
if ~isempty(excludecoordidx)
    excludecoord=varargin{excludecoordidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(stateconstraintidx)
    stateconstraintflag=varargin{stateconstraintidx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if isempty(stopcriterion)
    stopcriterion=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if ~isempty(excludecoordinate4epidx)
    excludecoordinate4ep=varargin{excludecoordinate4epidx+1};
end
if isempty(pathtype)
    pathtype='s';
end
targetvalue=targetvalue(:);

OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if exogenousfunction
    OCMATAE.exogenousinitialstates=exogenousinitialstates;
    OCMATAE.exogenousnumberofstates=length(exogenousinitialstates);
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianasymptotic=funch{6}{2};
OCMATAE.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATAE.hybridinfo=funch{7}{1};
OCMATAE.domain=funch{7}{2};
OCMATAE.guard=funch{7}{3};
OCMATAE.reset=funch{7}{4};
OCMATAE.switchtime=funch{7}{5};
OCMATAE.jacobianguard=funch{7}{7};
OCMATAE.jacobianreset=funch{7}{8};
OCMATAE.domaindiscretization=funch{7}{9};
OCMATAE.timesettransformation=funch{7}{10};
if objectivevaluecalc
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    if ~isautonomous(ocObj)
        OCMATAE.objectivefunctionderivativetime=funch{8}{4};
    end
end
OCMATAE.userbc=userbc;
if ~isempty(userbc)
    OCMATAE.userbcfunc=funch{5}{15};
end
if ~isautonomous(ocObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
end
if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamics=funch{4}{1};
    OCMATAE.exogenousjacobian=funch{4}{2};
    OCMATAE.exogenousparameterjacobian=funch{4}{3};
end

% function for workspace initialization
try
    OCMATAE.probleminit=funch{10};
catch
    OCMATAE.probleminit=[];
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

OCMATAE.stopcriterion=stopcriterion;

ocEParcarg=arcargument(ocEP);
ocEParcindex=arcarg2arcindex(ocEParcarg);

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specific model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

if isempty(asymptoticapproximation) && isempty(infinity)
    timesettransformation=OCMATAE.timesettransformation();
elseif ~isempty(infinity)
    timesettransformation.normalization=1;
    timesettransformation.infinity=infinity;
    timesettransformation.asymptoticapproximation=asymptoticapproximation;
else
    timesettransformation.normalization=1;
    timesettransformation.infinity=0;
    timesettransformation.asymptoticapproximation=asymptoticapproximation;
end
% initialize solution using the repeated entries of the equilibrium
sol=generateodestruct(ocEP,TrivialArcMeshNum,timesettransformation);
OCMATAE.freeendtime=freeendtime;
if freeendtime
    sol.parameters=[sol.parameters 0];
    OCMATAE.freeendtimecoord=length(sol.parameters);
end
numberofodes=canonicalsystemdimension(ocObj);

% add continuation parameter value
sol.parameters=[sol.parameters zeros(1,size(freevector,2)+1)];


if objectivevaluecalc
    o=objectivefunction(ocObj,sol);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
    numberofodes=numberofodes+1;
end
if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATAE.exogenousnumberofstates;
else
    OCMATAE.exogenousdynamicscoordinate=[];
end

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(1).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(1).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(1).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(1).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(1).odedim+(1:domaindata(1).aedim);
end
numode=domaindata(arcarg2arcindex(arcargument(ocEP))).odedim;

% test if algebraic equations have to be considered
numae=domaindata(arcarg2arcindex(arcargument(ocEP))).aedim;

% reduce Jacobian to ODE part
J=linearization(ocEP);
OCMATAE.linearization=J;
if ~isempty(fixedcoord)
    J(fixedcoord,:)=[];
    J(:,fixedcoord)=[];
end
if ~isempty(excludecoord)
    J(excludecoord,:)=[];
    J(:,excludecoord)=[];
end
OCMATAE.fixedcoord=fixedcoord;
OCMATAE.excludecoord=excludecoord;
OCMATAE.excludecoordinate4ep=[];

% mode and path specific variables
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=size(sol.y,1);
end
OCMATAE.inftimetransformation=timesettransformation.infinity;
OCMATAE.excludecoordinate4ep=excludecoordinate4ep;

switch pathtype
    case {'s','sc','cs','sts','ws'}
        OCMATAE.truncationtime=timesettransformation.asymptoticapproximation;
    case {'u','uc','cu','stu','wu','os'}
        OCMATAE.truncationtime=-timesettransformation.asymptoticapproximation;
        sol.arcinterval=-sol.arcinterval;
end
OCMATAE.saddlepoint=dependentvar(ocEP);
OCMATAE.pathtype=pathtype;
if isempty(asymptoticmatrix)
    asymptoticmatrix=asymptoticbc(J,pathtype,'c',ZeroDeviationTolerance);
end
if ~isempty(fixedcoord)
    OCMATAE.asymptoticmatrix=zeros(numode,size(asymptoticmatrix,2));
    coord=1:numode;
    coord(fixedcoord)=[];
    OCMATAE.asymptoticmatrix(coord,:)=asymptoticmatrix;
elseif ~isempty(excludecoord)
    OCMATAE.asymptoticmatrix=zeros(numode,size(asymptoticmatrix,2));
    coord=1:numode;
    coord(excludecoord)=[];
    OCMATAE.asymptoticmatrix(coord,:)=asymptoticmatrix;
else
    OCMATAE.asymptoticmatrix=asymptoticmatrix;
end
% if numae
%     % remove 
%     [abseigval,sortidx]=sort(abs(infoStruct.eigval));
%     asymptoticmatrix(:,infoStruct.index(sortidx(1:numae)))=[];
% end
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=OCMATAE.saddlepoint(contcoordinate);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
OCMATAE.implicitcontrolindex=domaindata(ocEParcindex).implicitcontrolindex;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATAE.freevector=freevector;
OCMATAE.freevectorindex=1:size(freevector,2);

OCMATAE.movinghorizon=0;

% test if pure state constraints are defined
if isempty(~stateconstraintflag)
    OCMATAE.stateconstraint=stateconstraint(ocObj);
else
    OCMATAE.stateconstraint=stateconstraintflag;
end
if OCMATAE.stateconstraint
    OCMATAE.jumpcostatecoord=[];
    OCMATAE.jumpcostateindex=[];
end
OCMATAE.objectivevaluecalc=objectivevaluecalc;

OCMATCONT.HE.numinitialcondition=numel(contcoordinate)+numae;
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATAE.ODEcoord=1:(size(sol.y,1)+OCMATAE.exogenousnumberofstates);
