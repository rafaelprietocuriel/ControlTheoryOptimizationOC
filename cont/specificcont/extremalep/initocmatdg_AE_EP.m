function sol=initocmatdg_AE_EP(dgObj,ocEP,contcoordinate,targetvalue,opt,varargin) 
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
sol=[];
pathtype='';
asymptoticapproximation=[];
freeendtime=[];
objectivevaluecalc=[];
asymptoticmatrix=[];
targetstate=[];
targetstatecoordinate=[];

if isempty(dgObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocEP)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if nargin==4
    opt=defaultocoptions;
end
asymptoticapproximationidx=find(strcmpi(varargin,'truncationtime'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'));
targetstateidx=find(strcmpi(varargin,'targetstate'));
targetstatecoordinateidx=find(strcmpi(varargin,'targetstatecoordinate'));

if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(asymptoticapproximationidx)
    asymptoticapproximation=varargin{asymptoticapproximationidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(pathtypeidx)
    pathtype=varargin{pathtypeidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(targetstateidx)
    targetstate=varargin{targetstateidx+1};
end
if ~isempty(targetstatecoordinateidx)
    targetstatecoordinate=varargin{targetstatecoordinateidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(pathtype)
    pathtype='s';
end
targetvalue=targetvalue(:);

OCMATAE.targetstate=targetstate;
OCMATAE.targetstatecoordinate=targetstatecoordinate;

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(dgObj);
OCMATCONT.modelfunc=modelspecificfunc(dgObj,'4SaddlePathContinuation');

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
OCMATAE.guard=funch{7}{3};
OCMATAE.reset=funch{7}{4};
if objectivevaluecalc
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    if ~isautonomous(dgObj)
        OCMATPS.objectivefunctionderivativetime=funch{8}{4};
    end
end
if ~isautonomous(dgObj)
    OCMATPS.canonicalsystemderivativetime=funch{2}{3};
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

% initialize solution using the repeated entries of the equilibrium
sol=generateodestruct(ocEP,TrivialArcMeshNum);
OCMATAE.freeendtime=freeendtime;
OCMATAE.playernum=playernum(dgObj);

if freeendtime
    sol.parameters=[sol.parameters 0];
    OCMATAE.freeendtimecoord=length(sol.parameters);
end
sol.parameters=[sol.parameters 0];
if objectivevaluecalc
    o=objectivefunction(dgObj,sol);
    OCMATAE.objectivevaluecoord=size(sol.y,1)+(1:OCMATAE.playernum);
    sol.y(end+(1:OCMATAE.playernum),:)=[zeros(OCMATAE.playernum,1) cumsum((o(:,1:end-1)+o(:,2:end))/2.*repmat(diff(time(dgObj,sol,1)),OCMATAE.playernum,1))];
end

% reduce Jacobian to ODE part
J=linearization(ocEP);
OCMATAE.linearization=J;

% mode and path specific variables
OCMATAE.parametervalue=parametervalue(dgObj);
OCMATAE.autonomous=isautonomous(dgObj);

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;

switch pathtype
    case {'s','sc','cs','sts','ws'}
        OCMATAE.truncationtime=asymptoticapproximation;
    case {'u','uc','cu','stu','wu','os'}
        OCMATAE.truncationtime=-asymptoticapproximation;
        sol.arcinterval=-sol.arcinterval;
end
OCMATAE.saddlepoint=dependentvar(ocEP);
OCMATAE.pathtype=pathtype;
if isempty(asymptoticmatrix)
    asymptoticmatrix=asymptoticbc(J,pathtype,'c',ZeroDeviationTolerance);
end
OCMATAE.asymptoticmatrix=asymptoticmatrix;
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=OCMATAE.saddlepoint(contcoordinate);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATAE.movinghorizon=0;

% test if pure state constraints are defined
OCMATAE.objectivevaluecalc=objectivevaluecalc;

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATCONT.numeq=size(sol.y,1);
OCMATCONT.equilibriumcoord=1:length(dependentvar(ocEP));
function sol=generateodestruct(ocEP,initnummesh)

sol=[];
if isempty(ocEP)
    return
end

sol.x=linspace(0,1,initnummesh);
sol.y=repmat(dependentvar(ocEP),1,initnummesh);
sol.parameters=[];
sol.arcinterval=[0 1]; %
sol.arcposition=[1;size(sol.y,2)]; %
sol.timehorizon=inf;
sol.arcarg=arcargument(ocEP);
sol.x0=sol.x(1);
sol.idata.tangent=[];
sol.idata.coeff=[];