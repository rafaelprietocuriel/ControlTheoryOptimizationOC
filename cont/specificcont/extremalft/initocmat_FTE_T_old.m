function sol=initocmat_FTE_T(ocObj,ocTrj,timeindex,targetvalue,initialcoordinate,varargin)
%
% initocmat_FTE_T initialization for asymptotic extremal calculation
%
% SOL=initocmat_FTE_T(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see initocmat_FTE_T) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
initialstate=[];
initialcostate=[];
initialendtime=[];
initialarcargument=[];
findoptimalhorizon=[];
hitstatevalue=[];
hitstatecoordinate=[];
fixendstate=[];
fixendcostate=[];
fixinitstate=[];
fixinitcostate=[];
objectivevaluecalc=[];
exogenousfunction=[];
jumpargument=[];
entryindex=[];
opt=[];
optimalhorizon=false;
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
optionidx=find(strcmpi(varargin,'option'));
initialstateidx=find(strcmpi(varargin,'initialstate'));
initialcostateidx=find(strcmpi(varargin,'initialcostate'));
initialendtimeidx=find(strcmpi(varargin,'initialendtime'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
fixinitcostateidx=find(strcmpi(varargin,'fixinitcostate'));
fixendcostateidx=find(strcmpi(varargin,'fixendcostate'));
optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
findoptimalhorizonidx=find(strcmpi(varargin,'findoptimalhorizon'));
hitstatevalueidx=find(strcmpi(varargin,'hitstatevalue'));
hitstatecoordinateidx=find(strcmpi(varargin,'hitstatecoordinate'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));

if ~isempty(fixendcostateidx)
    fixendcostate=varargin{fixendcostateidx+1};
end
if ~isempty(fixinitcostateidx)
    fixinitcostate=varargin{fixinitcostateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(jumpargumentidx)
    jumpargument=varargin{jumpargumentidx+1};
end
if ~isempty(entryindexidx)
    entryindex=varargin{entryindexidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
    initialcoordinate=[];
else
    if ~isempty(fixendstateidx)
    initialcoordinate=[];
    else
        fixinitstate=initialcoordinate;
    end
end
if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(findoptimalhorizonidx)
    findoptimalhorizon=varargin{findoptimalhorizonidx+1};
end
if ~isempty(hitstatevalueidx)
    hitstatevalue=varargin{hitstatevalueidx+1};
end
if ~isempty(hitstatecoordinateidx)
    hitstatecoordinate=varargin{hitstatecoordinateidx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
if isempty(findoptimalhorizon)
    findoptimalhorizon=false;
end
if isempty(hitstatevalue)
    hitstatevalue=[];
    hitstatecoordinate=[];
end

%%% generate initial octrajectory if not provided
if isempty(ocTrj)
    if ~isempty(initialstateidx)
        initialstate=varargin{initialstateidx+1};
        initialstate=initialstate(:);
    end
    if ~isempty(initialcostateidx)
        initialcostate=varargin{initialcostateidx+1};
        initialcostate=initialcostate(:);
    end
    if ~isempty(initialendtimeidx)
        initialendtime=varargin{initialendtimeidx+1};
    end
    if ~isempty(initialarcargumentidx)
        initialarcargument=varargin{initialarcargumentidx+1};
    end
    if ~isempty(optionidx)
        opt=varargin{optionidx+1};
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    if isempty(initialarcargument)
        initialarcargument=0;
    end
    if isempty(initialendtime)
        initialendtime=0;
    end
    n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocTrj.x=linspace(0,1,n);
    initocPt.y=initialstate;
    initocPt.x=0;
    initocPt.arcarg=initialarcargument;
    initocPt.arcinterval=[0 initialendtime];
    initocPt.arcposition=[1;1];
    if isempty(initialcostate)
        initialcostate=transversalitycondition(ocObj,octrajectory(initocPt));
    end
    ocTrj.y=[initialstate(:);initialcostate];
    ocTrj.y=ocTrj.y(:,ones(1,n));
    ocTrj.arcposition=[1;n];
    ocTrj.arcinterval=[0 initialendtime];
    ocTrj.arcarg=initialarcargument;
    ocTrj.x0=0;
    ocTrj.timehorizon=initialendtime;
    ocTrj.modelparameter=parametervalue(ocObj);
    ocTrj.modelname=modelname(ocObj);
    ocTrj=octrajectory(ocTrj);
end

OCMATFTE.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATFTE.exogenousinitialstates=exogenousinitialstates;
    OCMATFTE.exogenousnumberofstates=length(exogenousinitialstates);
end

% test if pure state constraints are defined
OCMATFTE.stateconstraint=stateconstraint(ocObj)&&~isempty(jumpargument);
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4FiniteHorizonPathContinuation');

% initialize global variable (OCMATFTE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATFTE.canonicalsystem=funch{1};
OCMATFTE.canonicalsystemjacobian=funch{2}{1};
OCMATFTE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATFTE.canonicalsystemhessian=funch{3}{1};
OCMATFTE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATFTE.bcinitial=funch{5}{1};
OCMATFTE.bctransversality=funch{5}{2};
OCMATFTE.bcoptimalhorizon=funch{5}{3};
if OCMATFTE.stateconstraint
    OCMATFTE.bcstateconstraint=funch{5}{7};
    OCMATFTE.bctransversalitysc=funch{5}{8};
end

% function for Jacobian
OCMATFTE.bcjacobianinitial=funch{6}{1};
OCMATFTE.bcjacobiantransversality=funch{6}{2};

% function describing the hybrid structure of the problem
OCMATFTE.hybridinfo=funch{7}{1};
OCMATFTE.domain=funch{7}{2};
OCMATFTE.guard=funch{7}{3};
OCMATFTE.reset=funch{7}{4};
OCMATFTE.switchtime=funch{7}{5};
OCMATFTE.jacobianguard=funch{7}{7};
OCMATFTE.jacobianreset=funch{7}{8};
OCMATFTE.domaindiscretization=funch{7}{9};

OCMATFTE.findoptimalhorizon=findoptimalhorizon;

if objectivevaluecalc
    OCMATFTE.objectivefunction=funch{8}{1};
    OCMATFTE.objectivefunctionjacobian=funch{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATFTE.objectivefunctionderivativetime=funch{8}{4};
    OCMATFTE.salvagevalue=funch{5}{6};
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamics=funch{4}{1};
    OCMATFTE.exogenousjacobian=funch{4}{2};
    OCMATFTE.exogenousparameterjacobian=funch{4}{3};
end

% general function
OCMATFTE.plotcontinuation=funch{11};
OCMATFTE.testadmissibility=funch{12};
OCMATFTE.datapath=funch{20};
OCMATFTE.saveintermediatefiles=funch{21};

hybridinfo=OCMATFTE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATFTE.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocTrj,solver(ocTrj));
OCMATFTE.transversalityconditioncs=0;

if OCMATFTE.stateconstraint
    OCMATFTE.jumpargument=jumpargument;
    OCMATFTE.entryindex=zeros(1,length(sol.arcinterval));
    OCMATFTE.entryindex(entryindex)=1:length(jumpargument);
    if ~isempty(jumpargument)
        entrytimecoordinate=length(sol.parameters)+1;
        sol.parameters=[sol.parameters jumpargument];
        OCMATFTE.entrytimecoordinate=entrytimecoordinate:length(sol.parameters);
    end
    if entryindex(end)==length(sol.arcinterval)
        OCMATFTE.transversalityconditioncs=1;
    end
end

arctimes=arcinterval(ocTrj);

numswitchtimes=length(sol.arcinterval(2:end));
switchtimecoord=length(sol.parameters)+1;
if isempty(sol.parameters)
    sol.parameters=[sol.parameters sol.arcinterval(2:end-1)];
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=[sol.parameters sol.arcinterval(2:end-1)];
end

sol.parameters=[sol.parameters arctimes(timeindex)];
OCMATFTE.switchtimecoord=switchtimecoord:length(sol.parameters);
if optimalhorizon
    sol.parameters=[sol.parameters 0];
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
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj) || objectivevaluecalc>1
    OT=discountedsalvagevalue(ocObj,ocTrj);
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(end+1,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
end

% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);

OCMATFTE.initialtime=sol.x0;
%OCMATFTE.initialstate=depvar(initialcoordinate,1);
OCMATFTE.initialcoordinate=initialcoordinate;
OCMATFTE.varyparameterindex=timeindex;
OCMATFTE.targetvalue=targetvalue;
OCMATFTE.optimalhorizon=optimalhorizon;

if optimalhorizon
    hamiltonian(ocObj,ocTrj,1);
    OCMATFTE.startvalue=ans(end);
    OCMATFTE.continuationvector=-OCMATFTE.startvalue;
end
OCMATFTE.objectivevaluecalc=objectivevaluecalc;
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;
OCMATFTE.objectivevaluecoord=[];
if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATFTE.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.exogenousnumberofstates;
else
    OCMATFTE.exogenousdynamicscoordinate=[];
end

OCMATFTE.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATFTE.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATFTE.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATFTE.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATFTE.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATFTE.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
if objectivevaluecalc
    OCMATFTE.dFDO=zeros(dimensioncanonicalsystem,1);
    OCMATFTE.dFODO=0;
    OCMATFTE.dFODX=zeros(1,dimensioncanonicalsystem);
    if exogenousfunction
        OCMATFTE.dFODE=zeros(1,OCMATFTE.exogenousnumberofstates);
    end
    OCMATFTE.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATFTE.exogenousfunction
    OCMATFTE.dFDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDE=zeros(OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDX=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    if OCMATFTE.objectivevaluecalc
        OCMATFTE.dFEDO=zeros(OCMATFTE.exogenousnumberofstates,1);
    end
    OCMATFTE.dFEDPAR=zeros(OCMATFTE.exogenousnumberofstates,length(sol.parameters));
end
OCMATFTE.Jpar=[OCMATFTE.dFDPAR;OCMATFTE.dFODPAR;OCMATFTE.dFEDPAR];
OCMATFTE.ODEcoord=1:size(OCMATFTE.Jpar,1);
 
OCMATFTE.Jext=[[OCMATFTE.dFDO OCMATFTE.dFDE];[OCMATFTE.dFODO OCMATFTE.dFODE];[OCMATFTE.dFEDO OCMATFTE.dFEDE]];

OCMATFTE.statenum=statenum(ocObj);
OCMATFTE.fixinitstatecoord=fixinitstate;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitcostatecoord=fixinitcostate;
OCMATFTE.fixendcostatecoord=fixendcostate;
OCMATFTE.initstate=depvar(fixinitstate,1);
OCMATFTE.endstate=depvar(fixendstate,end);
OCMATFTE.endcostate=depvar(OCMATFTE.statenum+fixendstate,end);
OCMATFTE.hitstatevalue=hitstatevalue;
OCMATFTE.hitstatecoordinate=hitstatecoordinate;

%OCMATFTE.continitstate=false;
% if isempty(fixinitstate)
%     OCMATFTE.startvalue=depvar(fixinitstate,1);
%     %OCMATFTE.continitstate=true;
% elseif isempty(fixendstateidx)
%     OCMATFTE.startvalue=depvar(fixendstate,end);
% end
% if ~isempty(fixendstate)
%     OCMATFTE.endstate=depvar(fixendstate,end);
%     OCMATFTE.initstate=[];
% elseif ~isempty(fixinitstate)
%     OCMATFTE.initstate=depvar(fixinitstate,1);
%     OCMATFTE.endstate=[];
% end
pathname=OCMATFTE.datapath();
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
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