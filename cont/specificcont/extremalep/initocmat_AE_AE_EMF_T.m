function sol=initocmat_AE_AE_EMF_T(ocObj,ocAsym,timeindex,initialcoordinate,opt,varargin)
%
% INITOCMAT_AE_AE initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_AE(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_AE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
jump=[]; % provide possible jumps of costate for pure state constraints
pathtpe=pathtype(ocAsym);
explicitemfcoord=[];
independentemfcoord=[];
constjac=[];
targetvalue=[];
targettype='';
if isempty(pathtpe)
    pathtpe='s';
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==5
    opt=defaultocoptions;
end
constjacidx=find(strcmpi(varargin,'constantjacobian'),1);
explicitemfcoordidx=find(strcmpi(varargin,'explicitcoordinate'),1);
indemfcoordidx=find(strcmpi(varargin,'indepcoordinate'),1);
targettypeidx=find(strcmpi(varargin,'targettype'));
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if isempty(targetvalue)
    targetvalue=1e-1;
end
if isempty(targettype)
    targettype='T';
end
if ~isempty(explicitemfcoordidx)
    explicitemfcoord=varargin{explicitemfcoordidx+1};
end
if ~isempty(constjacidx)
    constjac=varargin{constjacidx+1};
end
if ~isempty(indemfcoordidx)
    independentemfcoord=varargin{indemfcoordidx+1};
end

solverInfoStruct=solverinfo(ocAsym);
if isfield(solverInfoStruct,'explicitemfcoord') && isempty(explicitemfcoordidx)
    explicitemfcoord=solverInfoStruct.explicitemfcoord;
end
if isfield(solverInfoStruct,'independentemfcoord') && isempty(indemfcoordidx)
    independentemfcoord=solverInfoStruct.independentemfcoord;
end
if isfield(solverInfoStruct,'constantjacobian') && isempty(constjacidx)
    constjac=solverInfoStruct.constantjacobian;
end
if isempty(constjac)
    constjac=0;
end
targetvalue=targetvalue(:);
stableflag=~isempty(strfind(pathtpe,'s'));

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
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

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocAsym,solver(ocAsym));
arctimes=arcinterval(ocAsym);
if isempty(timeindex)
    timeindex=length(arctimes);
end

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
numswitchtimes=length(sol.arcinterval(2:end-1));
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% mode and path specific variables
limSet=limitset(ocAsym);
JacobianMatrix=linearization(limSet);
limSetdepvar=dependentvar(limSet);
depvar=dependentvar(ocAsym);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
[OCMATAE.asymptoticmatrix, OCMATAE.numstable OCMATAE.numunstable OCMATAE.numcenter]=asymptoticbc(JacobianMatrix,pathtpe,'c',ZeroDeviationTolerance);

OCMATAE.initialtime=sol.x0;
OCMATAE.inftimetransformation=inftimetransformation(ocAsym);
if isempty(OCMATAE.inftimetransformation)
    OCMATAE.inftimetransformation=0;
end
switch pathtpe
    case 's'
        OCMATAE.truncationtime=sol.arcinterval(end);
        if ~constjac
            Y=zeros(OCMATAE.numunstable+OCMATAE.numcenter,OCMATAE.numstable);
            OCMATAE.Y=Y;
            OCMATAE.Q0=computeBase(JacobianMatrix,stableflag,OCMATAE.numstable);
            dimSubSpace=OCMATAE.numstable;
            OCMATCONT.HE.Ycoord=reshape(length(limSetdepvar)+(1:OCMATAE.numstable*(OCMATAE.numunstable+OCMATAE.numcenter)),OCMATAE.numstable,OCMATAE.numunstable+OCMATAE.numcenter);
            OCMATAE.Id=eye(OCMATAE.numunstable+OCMATAE.numcenter);
            OCMATAE.dimSubSpace=dimSubSpace;
            OCMATAE.subspacedim=OCMATAE.dimSubSpace;
            OCMATAE.numY=numel(Y);
            OCMATAE.orthspacedim=length(limSetdepvar)-OCMATAE.subspacedim;
        end
    case 'u'
        OCMATAE.truncationtime=-sol.arcinterval(end);
end
OCMATAE.linearization=JacobianMatrix;
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.explicitemfcoord=explicitemfcoord; % a user defined function has to provide the explicit values manifold equilibrium
OCMATAE.independentemfcoord=independentemfcoord; % coordinates used as independent variables for the equilibrium manifold
OCMATAE.dependentemfcoord=setdiff(1:length(limSetdepvar),[explicitemfcoord independentemfcoord]); % number of dependent manifold equations.
OCMATAE.emfcoord=setdiff(1:length(limSetdepvar),explicitemfcoord);
OCMATAE.emfindex=1:length(OCMATAE.emfcoord);
OCMATAE.initialstate=depvar(initialcoordinate,1);
OCMATAE.initialcoordinate=initialcoordinate;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=false;
OCMATAE.constantjacobian=constjac;

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATAE.stableflag=stableflag;
OCMATAE.targetvalue=targetvalue;
OCMATAE.targettype=targettype;

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;

if isempty(sol.parameters)
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1)];
elseif length(sol.parameters)==numswitchtimes
    sol.parameters=[jump(:).' sol.parameters];
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1)];
end
if ~constjac
    sol.parameters=[limSetdepvar(OCMATAE.emfcoord).' Y(:).' sol.parameters];
    OCMATAE.switchtimecoord=length(limSetdepvar)+length(Y(:))+(1:numswitchtimes);
else
    sol.parameters=[limSetdepvar(OCMATAE.emfcoord).' sol.parameters];
    OCMATAE.switchtimecoord=length(limSetdepvar)+(1:numswitchtimes);
end
sol.parameters=[sol.parameters arctimes(timeindex)];

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=parameters(ocTrj);
if isfield(ocTrj.solverinfo,'conttype')
    switch ocTrj.solverinfo.conttype
        case 'extremal2ep'
            if ~isempty(sol.parameters)
                numcontpar=length(continuationparameter(ocTrj));
                sol.parameters(end-numcontpar+1:end)=[];
            end
        otherwise
            sol.parameters=[];
    end
end
sol.solver=solvername;

if isfield(ocTrj.solverinfo,'tangent')
    sol.solverinfo.tangent=ocTrj.solverinfo.tangent;
else
    sol.solverinfo.tangent=[];
end
if isfield(ocTrj.solverinfo,'coeff')
    sol.solverinfo.coeff=ocTrj.solverinfo.coeff;
else
    sol.solverinfo.coeff=[];
end
if isfield(ocTrj.solverinfo,'yp')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ypmid')
    sol.solverinfo.ypmid=ocTrj.solverinfo.ypmid;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.ymid=ocTrj.solverinfo.ymid;
end