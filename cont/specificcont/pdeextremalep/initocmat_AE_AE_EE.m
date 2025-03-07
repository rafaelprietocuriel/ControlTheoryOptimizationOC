function sol=initocmat_AE_AE_EE(ppdeObj,pdeAsym,contcoordinate,targetvalue,opt,varargin) 
%
% initocmat_AE_AE_EE initialization for asymptotic extremal calculation
% converging to a solution of an elliptic PDE.
%
% SOL=initocmat_AE_AE_EE(PPDEOBJ,PDEASYM,CONTCOORDINATE,TARGETVALUE) a stable
% saddle path calculation is initialized. 
% PPDEOBJ          ... corresponding optimal control model
% PDEASYM          ... asymptotic extremal converging to an elliptic
%                      equilibrium
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
% initial solution is given as well. In that sense the OCEE argument
% performs two tasks. The searched for solution converges to OCEE and OCEE
% is the initial solution (mu=0) of the continuation process.
%
% The denomination as TARGETVALUE maybe misleading but from the
% continuation point of view it denotes the target. From the problem
% perspective it denotes the initial state values of the searched for
% solution. 
% 
% SOL=INITOCMAT_AE_EP(PPDEOBJ,OCEE,CONTCOORDINATE,TARGETVALUE,OPT) with the
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
sol=[];
pathtpe='';
movinghorizon=[];
freevector=[];
freeendtime=[];
fixedcoord=[];
excludecoord=[];
objectivevaluecalc=[];
asymptoticmatrix=[];
pathtpe=pathtype(pdeAsym);

if isempty(ppdeObj)
    ocmatmsg('ppde model is empty.\n')
    return
end

pdeEP=limitset(pdeAsym);
if isempty(pdeEP)
    ocmatmsg('pde equilibrium is empty.\n')
    return
end
if nargin==4
    opt=defaultocoptions;
end
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
freevectoridx=find(strcmpi(varargin,'freevector'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
excludecoordidx=find(strcmpi(varargin,'excludecoord'),1);
fixedcoordidx=find(strcmpi(varargin,'fixedcoord'),1);
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'));

if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(pathtypeidx)
    pathtpe=varargin{pathtypeidx+1};
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
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(movinghorizon)
    movinghorizon=0;
end
if isempty(pathtpe)
    pathtpe='s';
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ppdeObj);
OCMATCONT.modelfunc=modelspecificfunc(ppdeObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
if ~isautonomous(ppdeObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
end
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
    if ~isautonomous(ppdeObj)
        OCMATAE.objectivefunctionderivativetime=funch{8}{4};
    end
end

OCMATAE.autonomous=isautonomous(ppdeObj);

if ~OCMATAE.autonomous
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
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

OCMATAE.saddlepoint=dependentvar(pdeEP);
OCMATAE.femdata=femdata(pdeAsym);
% initialize solution using the repeated entries of the equilibrium
sol=generatesolstruct(pdeAsym);
OCMATAE.freeendtime=freeendtime;
if freeendtime
    sol.parameters=[sol.parameters 0];
    OCMATAE.freeendtimecoord=length(sol.parameters);
end

% add continuation parameter value
sol.parameters=[sol.parameters zeros(1,size(freevector,2)+1)];


if objectivevaluecalc
    sol.femdata=OCMATAE.femdata;
    pdeTrj=pdetrajectory(sol);
    OCMATAE.objectivevaluecoord=size(sol.y,1)+1;
    o=objectivefunction(ppdeObj,pdeTrj);
    sol.y(end+(1:OCMATAE.femdata.gridnum),:)=[zeros(OCMATAE.femdata.gridnum,1) cumsum((o(:,1:end-1)+o(:,2:end)))/2.*repmat(diff(sol.x*OCMATAE.truncationtime),OCMATAE.femdata.gridnum,1)];
    OCMATAE.objectivevaluecoord=OCMATAE.objectivevaluecoord:size(sol.y,1);
end

pdeAsymarcarg=arcargument(pdeAsym);
for ii=1:numel(pdeAsymarcarg)
    OCMATCONT.DOMAINDDATA(ii).numode=length(OCMATAE.saddlepoint);
    OCMATCONT.DOMAINDDATA(ii).numae=0;
    OCMATCONT.DOMAINDDATA(ii).daeorder=[];
    OCMATCONT.DOMAINDDATA(ii).numeq=length(OCMATAE.saddlepoint);
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=OCMATCONT.DOMAINDDATA(ii).numode+OCMATAE.femdata.gridnum;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(1).numeq+OCMATAE.femdata.gridnum;%number of equations
    end
    numode=OCMATCONT.DOMAINDDATA(ii).numode;
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(1).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:OCMATCONT.DOMAINDDATA(ii).numode;
    OCMATCONT.DOMAINDDATA(ii).aecoord=[];
end

% reduce Jacobian to ODE part
J=full(linearization(pdeEP));
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
OCMATAE.parametervalue=parametervalue(ppdeObj);

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;
switch pathtpe
    case {'s','sc','cs','sts','ws'}
        OCMATAE.truncationtime=sol.arcinterval(end);
    case {'u','uc','cu'}
        OCMATAE.truncationtime=-abs(sol.arcinterval(end));
        sol.arcinterval=-abs(sol.arcinterval(end));
end
OCMATAE.pathtype=pathtpe;
if isempty(asymptoticmatrix)
    asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
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

depVar=dependentvar(pdeAsym);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depVar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATAE.freevector=freevector;
OCMATAE.freevectorindex=1:size(freevector,2);

OCMATAE.movinghorizon=0;
OCMATAE.objectivevaluecalc=objectivevaluecalc;

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;


function sol=generatesolstruct(pdeAsym)

sol.x=independentvar(pdeAsym);
sol.y=dependentvar(pdeAsym);
sol.arcarg=arcargument(pdeAsym);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(pdeAsym);
sol.arcposition=arcposition(pdeAsym);
sol.parameters=[];
