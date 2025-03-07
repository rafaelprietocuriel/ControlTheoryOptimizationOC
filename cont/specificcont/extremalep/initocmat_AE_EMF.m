function sol=initocmat_AE_EMF(ocObj,ocEP,contcoordinate,targetvalue,opt,varargin) 
%
% INITOCMAT_AE_EMF initialization for asymptotic extremal calculation
% converging to an equilibrium manifold
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE

sol=[];
pathtpe='';
asymptoticapproximation=[];
infinity=[];
explicitemfcoord=[];
independentemfcoord=[];
exogenousfunction=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocEP)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if ~testconsistency(ocEP,ocObj)
    ocmaterror('The equilibrium and oc model are not consistent.')
end
if nargin==4
    opt=defaultocoptions;
end
asymptoticapproximationidx=find(strcmpi(varargin,'truncationtime'));
infinityidx=find(strcmpi(varargin,'inftimetransformation'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
constjac=find(strcmpi(varargin,'constantjacobian'));
explicitemfcoordidx=find(strcmpi(varargin,'explicitcoordinate'),1);
indemfcoordidx=find(strcmpi(varargin,'indepcoordinate'),1);
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));

if ~isempty(asymptoticapproximationidx)
    asymptoticapproximation=varargin{asymptoticapproximationidx+1};
end
if ~isempty(infinityidx)
    infinity=varargin{infinityidx+1};
end
if ~isempty(infinity) && infinity
    asymptoticapproximation=1;
end
if ~isempty(pathtypeidx)
    pathtpe=varargin{pathtypeidx+1};
end
if ~isempty(indemfcoordidx)
    independentemfcoord=varargin{indemfcoordidx+1};
end
if ~isempty(explicitemfcoordidx)
    explicitemfcoord=varargin{explicitemfcoordidx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end

if isempty(pathtpe)
    pathtpe='s';
end
if isempty(constjac)
    constjac=0;
end

OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if exogenousfunction
    OCMATAE.exogenousinitialstates=exogenousinitialstates;
    OCMATAE.exogenousnumberofstates=length(exogenousinitialstates);
end

targetvalue=targetvalue(:);
stableflag=~isempty(strfind(pathtpe,'s'));

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

if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamics=funch{4}{1};
    OCMATAE.exogenousjacobian=funch{4}{2};
    OCMATAE.exogenousparameterjacobian=funch{4}{3};
end

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.equilibrium=funch{5}{4};
OCMATAE.explicitequilibriumvalue=funch{5}{5};

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

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

ocEParcarg=arcargument(ocEP);
ocEParcindex=arcarg2arcindex(ocEParcarg);

OCMATAE.movinghorizon=0;

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

% add continuation parameter value
sol.parameters=0;
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
%numode=domaindata(arcarg2arcindex(arcargument(ocEP))).odedim;

% test if algebraic equations have to be considered
numae=domaindata(arcarg2arcindex(arcargument(ocEP))).aedim;

limSetdepvar=dependentvar(ocEP);
JacobianMatrix=linearization(ocEP);
limSetDim=length(limSetdepvar);
% mode and path specific variables
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=[];
OCMATAE.saddlepoint=limSetdepvar;
OCMATAE.linearization=JacobianMatrix;
OCMATAE.pathtype=pathtpe;
OCMATAE.constantjacobian=constjac;
[asymptoticmatrix, OCMATAE.numstable OCMATAE.numunstable OCMATAE.numcenter infoStruct]=asymptoticbc(JacobianMatrix,pathtpe,'c',ZeroDeviationTolerance);
if isempty(independentemfcoord)
    independentemfcoord=setdiff(1:limSetDim,explicitemfcoord);
end

if OCMATAE.numcenter>length(independentemfcoord)
    ocmatmsg('The dimension of the null space exceeds the number of independent coordinates of the equilibrium manifold.')
end
if OCMATAE.numcenter<length(independentemfcoord)
    ocmatmsg('The number of independent coordinates of the equilibrium manifold exceeds the dimsion of the null space.')
end
OCMATAE.inftimetransformation=timesettransformation.infinity;
switch pathtpe
    case 's'
        OCMATAE.truncationtime=timesettransformation.asymptoticapproximation;
        if ~constjac
            OCMATAE.subspacedim=OCMATAE.numstable;
        end
    case 'sc'
        OCMATAE.truncationtime=timesettransformation.asymptoticapproximation;
        if ~constjac
            OCMATAE.subspacedim=OCMATAE.numstable+OCMATAE.numcenter;
        end
    case 'u'
        OCMATAE.truncationtime=-timesettransformation.asymptoticapproximation;
        if ~constjac
            OCMATAE.subspacedim=OCMATAE.numunstable;
        end
end
if ~constjac
    OCMATAE.orthspacedim=limSetDim-OCMATAE.subspacedim;
    Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);
    OCMATAE.Y=Y;
    OCMATAE.Q0=infoStruct.Q;
%         Q=computeBase(JacobianMatrix,~stableflag,OCMATAE.subspacedim);
    OCMATCONT.HE.Ycoord=reshape(limSetDim-length(explicitemfcoord)+(1:OCMATAE.orthspacedim*OCMATAE.subspacedim),OCMATAE.orthspacedim,OCMATAE.subspacedim);
    OCMATAE.Id=eye(OCMATAE.orthspacedim);
    OCMATAE.numY=numel(Y);
end
OCMATAE.explicitemfcoord=explicitemfcoord; % a user defined function has to provide the explicit values manifold equilibrium
OCMATAE.independentemfcoord=independentemfcoord; % coordinates used as independent variables for the equilibrium manifold
OCMATAE.dependentemfcoord=setdiff(1:limSetDim,[explicitemfcoord independentemfcoord]); % number of dependent manifold equations.
OCMATAE.emfcoord=setdiff(1:limSetDim,explicitemfcoord);
OCMATAE.emfindex=1:length(OCMATAE.emfcoord);
OCMATAE.asymptoticmatrix=asymptoticmatrix;
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=OCMATAE.saddlepoint(contcoordinate);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
OCMATAE.implicitcontrolindex=domaindata(ocEParcindex).implicitcontrolindex;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;

if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATAE.exogenousnumberofstates;
else
    OCMATAE.exogenousdynamicscoordinate=[];
end
% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
if OCMATAE.stateconstraint
    OCMATAE.jumpcostatecoord=[];
    OCMATAE.jumpcostateindex=[];
end
if ~constjac
    sol.parameters=[limSetdepvar(OCMATAE.emfcoord).' Y(:).' 0];
else
    sol.parameters=[limSetdepvar(OCMATAE.emfcoord).' 0];
end
% 
% if ~constjac
%     OCMATCONT.HE.Ycoord=OCMATCONT.HE.Ycoord;
% end
OCMATCONT.HE.numinitialcondition=numel(contcoordinate)+numae;
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);
OCMATAE.stableflag=stableflag;

OCMATCONT.codimension=1;

OCMATAE.ODEcoord=1:dimensioncanonicalsystem+OCMATAE.exogenousnumberofstates;
OCMATAE.freeparameter=[];