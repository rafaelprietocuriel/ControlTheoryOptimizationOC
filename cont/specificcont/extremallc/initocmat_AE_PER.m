function sol=initocmat_AE_PER(ocObj,ocPer,contcoordinate,targetvalue,opt,varargin)
%
% initocmat_AE_PER initialization for asymptotic extremal calculation
%
% SOL=initocmat_AE_PER(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
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

clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtype='';
asymptoticapproximation=[];
objectivevaluecalc=[];
fixedcoord=[];
asymptoticmatrix=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocPer)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if nargin==4
    opt=defaultocoptions;
end
asymptoticapproximationidx=find(strcmpi(varargin,'truncationtime'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
fixedcoordidx=find(strcmpi(varargin,'fixedcoord'),1);
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'),1);
if ~isempty(asymptoticapproximationidx)
    asymptoticapproximation=varargin{asymptoticapproximationidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(pathtypeidx)
    pathtype=varargin{pathtypeidx+1};
end
if isempty(pathtype)
    pathtype='s';
end
if ~isempty(fixedcoordidx)
    fixedcoord=varargin{fixedcoordidx+1};
end
if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
asymptoticbcmethod=getocoptions(opt,'OCCONTARG','AsymptoticBCMethod');
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
if ~isautonomous(ocObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATAE.objectivefunctionderivativetime=funch{8}{4};
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

OCMATAE.freeendtime=[];

ocPerarcarg=arcargument(ocPer);

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specific model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
numode=domaindata(arcarg2arcindex(arcargument(ocPer))).odedim;

asymptoticapproximation=ceil(asymptoticapproximation/period(ocPer));
timesettransformation.normalization=1;
timesettransformation.infinity=0;
timesettransformation.asymptoticapproximation=asymptoticapproximation;
% initialize solution using the repeated entries of the equilibrium
ocPer=removemonodromydynamics(ocPer);
sol=generateodestruct(ocPer,timesettransformation);
% add continuation parameter value
%sol.parameters=[sol.parameters 0];
if objectivevaluecalc
    sol.y(end+1,:)=objectivevalue(ocObj,sol);
end


% mode and path specific variables
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.fixedcoord=fixedcoord;
OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=size(sol.y,1);
end
OCMATAE.inftimetransformation=timesettransformation.infinity;
OCMATAE.truncationtime=timesettransformation.asymptoticapproximation*period(ocPer)+sol.x0;
OCMATAE.saddlepoint=dependentvar(ocPer);
OCMATAE.saddlepoint=OCMATAE.saddlepoint(:,end);
%OCMATAE.linearization=J;
OCMATAE.pathtype=pathtype;
OCMATAE.objectivevaluecalc=objectivevaluecalc;
OCMATAE.limitset=ocPer;

if isempty(asymptoticmatrix)
    % reduce Jacobian to ODE part
    J=jacobian(ocPer);
    J=J(1:numode,1:numode);
    if ~isempty(fixedcoord)
        J(fixedcoord,:)=[];
        J(:,fixedcoord)=[];
    end
    asymptoticmatrix=asymptoticbc(J,pathtype,'d',ZeroDeviationTolerance,asymptoticbcmethod);
end
if ~isempty(fixedcoord)
    OCMATAE.asymptoticmatrix=zeros(numode,size(asymptoticmatrix,2));
    coord=1:numode;
    coord(fixedcoord)=[];
    OCMATAE.asymptoticmatrix(coord,:)=asymptoticmatrix;
else
    OCMATAE.asymptoticmatrix=asymptoticmatrix;
end
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=OCMATAE.saddlepoint(contcoordinate);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
OCMATAE.implicitcontrolindex=[];%domaindata(ocPerarcindex).implicitcontrolindex;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

sol.parameters=[sol.parameters OCMATAE.truncationtime];
OCMATAE.truncationtimecoord=length(sol.parameters);
OCMATAE.freeendtime=true;
sol.parameters=[sol.parameters 0];

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATAE.stateconstraint=0;

function sol=generateodestruct(ocPer,timesettransformation)

sol=[];
if isempty(ocPer)
    return
end
if length(arcargument(ocPer))==1
    repmat(0:timesettransformation.asymptoticapproximation-1,length(independentvar(ocPer)),1);
    sol.x=repmat(independentvar(ocPer),1,timesettransformation.asymptoticapproximation)+ans(:)';
    sol.y=repmat(dependentvar(ocPer),1,timesettransformation.asymptoticapproximation);
    removeidx=find(diff(sol.x)==0);
    sol.x(removeidx)=[];
    sol.x=sol.x/sol.x(end);
    sol.y(:,removeidx)=[];
    sol.parameters=[];
    sol.arcinterval=[0 timesettransformation.asymptoticapproximation*period(ocPer)]; %
    sol.arcarg=arcargument(ocPer);
    sol.timehorizon=inf;
else
    t0=independentvar(ocPer);
    arcintv0=arcinterval(ocPer);
    y0=dependentvar(ocPer);
    sol.x=t0;
    sol.arcinterval=arcintv0;
    for ii=2:timesettransformation.asymptoticapproximation
        sol.x=[sol.x t0+sol.x(end)];
        sol.arcinterval=[sol.arcinterval arcintv0(2:end)+sol.arcinterval(end)];
    end
    sol.y=repmat(y0,1,timesettransformation.asymptoticapproximation);
    sol.arcarg=repmat(arcargument(ocPer),1,timesettransformation.asymptoticapproximation);
    arcposition=find(diff(sol.x)==0);
    sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
    sol=mergearc(sol);
end
sol.x0=initialtime(ocPer);
sol.arcinterval=sol.x0+sol.arcinterval;
arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.parameters=sol.arcinterval(2:end-1);
sol.idata.tangent=[];
sol.idata.coeff=[];