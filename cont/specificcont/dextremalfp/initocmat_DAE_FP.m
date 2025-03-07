function sol=initocmat_DAE_FP(docObj,ocFP,contcoordinate,targetvalue,opt,varargin) 
%
% INITOCMAT_DAE_FP initialization for asymptotic dextremal calculation
%
% SOL=INITOCMAT_DAE_FP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
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
% SOL=INITOCMAT_DAE_FP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE,OPT) with the
% option structure OPT the threshold 'ZeroDeviationTolerance' and initial
% number of discretization points 'TrivialArcMeshNum' for the equilibrium
% solution can be changed. 
%   ZeroDeviationTolerance ... provides the tolerance to classify an 
%                              eigenvalue numerically as zero.
%   TrivialArcMeshNum      ... provides the number of points for the
%                              constant solution at the equilbrium. 
%
% SOL=INITOCMAT_DAE_FP(...,'TruncationTime',T) the truncation of the
% infinite time horizon to the finite time T 
% 
% SOL=INITOCMAT_DAE_FP(...,'PathType',p)
%   p='s' (default) stable saddle-path calculation
%   p='u' unstable saddle-path calculation

clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtype='';
asymptoticapproximation=[];
if isempty(docObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocFP)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if ~testconsistency(ocFP,docObj)
    ocmaterror('The equilibrium and oc model are not consistent.')
end
if nargin==4
    opt=defaultocoptions;
end
asymptoticapproximationidx=find(strcmpi(varargin,'truncationtime'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
if ~isempty(asymptoticapproximationidx)
    asymptoticapproximation=varargin{asymptoticapproximationidx+1};
end
if isempty(asymptoticapproximation)
    asymptoticapproximation=100;
end
if ~isempty(pathtypeidx)
    pathtype=varargin{pathtypeidx+1};
end
if isempty(pathtype)
    pathtype='s';
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(docObj);
OCMATCONT.modelfunc=modelspecificfunc(docObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystemmap=funch{1};
OCMATAE.canonicalsystemmapjacobian=funch{2}{1};
OCMATAE.canonicalsystemmapparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemmaphessian=funch{3}{1};
OCMATAE.canonicalsystemmapparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianasymptotic=funch{6}{2};
OCMATAE.bcjacobiantransversality=funch{6}{3};

% general function
OCMATAE.findarcposition=funch{10};
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

% initialize solution using the repeated entries of the equilibrium
sol=generateodiffestruct(ocFP,asymptoticapproximation,pathtype);

% reduce Jacobian to ODE part
J=explicitjacobian(ocFP);

% mode and path specific variables
OCMATAE.parametervalue=parametervalue(docObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.truncationtime=asymptoticapproximation;
OCMATAE.saddlepoint=dependentvar(ocFP);
OCMATAE.saddlepoint=OCMATAE.saddlepoint(:,1);
OCMATAE.linearization=linearization(ocFP);
OCMATAE.pathtype=pathtype;
asymptoticmatrix=asymptoticbc(J,pathtype,'d',ZeroDeviationTolerance);

OCMATAE.asymptoticmatrix=asymptoticmatrix;
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=OCMATAE.saddlepoint(contcoordinate);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);
OCMATCONT.HE.numdvariables=numel(sol.y);
OCMATCONT.codimension=1;
OCMATCONT.numarc=1;

function sol=generateodiffestruct(ocFP,asymptoticapproximation,pathtype)

switch pathtype
    case 's'
        sol.x=0:asymptoticapproximation;
        sol.y=[ocFP.y0 repmat(ocFP.y,1,asymptoticapproximation)];
    case 'u'
        sol.x=asymptoticapproximation:-1:0;
        sol.y=[repmat(ocFP.y,1,asymptoticapproximation) ocFP.y0];
end
sol.x0=0;
% add continuation parameter value
sol.parameters=0;
sol.arcarg=arcargument(ocFP);
sol.arcposition=[1;size(sol.y,2)];