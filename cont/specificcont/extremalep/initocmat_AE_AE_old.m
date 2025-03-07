 function sol=initocmat_AE_AE(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
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
movinghorizon=[];
pathtpe=pathtype(ocAsym);
freevector=[];
excludecoord=[];
fixedcoord=[];
stopcriterion=[];
objectivevaluecalc=[];
asymptoticmatrix=[];
excludecoordinate4ep=[];
exogenousfunction=[];
userbc=[];

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
if nargin==4
    opt=defaultocoptions;
end
jumpidx=find(strcmpi(varargin,'jump'));
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
freevectoridx=find(strcmpi(varargin,'freevector'));
excludecoordidx=find(strcmpi(varargin,'excludecoord'),1);
fixedcoordidx=find(strcmpi(varargin,'fixedcoord'),1);
stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'));
excludecoordinate4epidx=find(strcmpi(varargin,'excludecoordinate4ep'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
userbcidx=find(strcmpi(varargin,'userbc'));

if ~isempty(excludecoordinate4epidx)
    excludecoordinate4ep=varargin{excludecoordinate4epidx+1};
end
if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(jumpidx)
    jump=varargin{jumpidx+1};
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if isempty(movinghorizon)
    movinghorizon=0;
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(excludecoordidx)
    excludecoord=varargin{excludecoordidx+1};
end
if ~isempty(fixedcoordidx)
    fixedcoord=varargin{fixedcoordidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if isempty(stopcriterion)
    stopcriterion=0;
end
if ~isempty(pathtypeidx)
    pathtpe=varargin{pathtypeidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end

targetvalue=targetvalue(:);

OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if exogenousfunction
    OCMATAE.exogenousinitialstates=exogenousinitialstates;
    OCMATAE.exogenousnumberofstates=length(exogenousinitialstates);
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
OCMATAE.freeendtime=[];

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
if ~isautonomous(ocObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
end
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
    OCMATAE.objectivefunctionderivativetime=funch{8}{4};
end
OCMATAE.userbc=userbc;
if ~isempty(userbc)
    OCMATAE.userbcfunc=funch{5}{15};
end

% if stopcriterion
%     OCMATAE.stopcriterionfunc=funch{5}{6};
% end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

OCMATAE.stopcriterion=stopcriterion;

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocAsym,solver(ocAsym));

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
sccounter=0;
if OCMATAE.stateconstraint && ~isempty(jump)
    arcarg4sc=arcargwithactivestateconstraint(ocObj);
    for ii=2:length(sol.arcarg)
        if ismember(sol.arcarg(ii),arcarg4sc)
            sccounter=sccounter+1;
            OCMATAE.jumpcostateindex(sccounter)=ii;
        end
    end
    if sccounter
        OCMATAE.jumpcostatecoord=1:sccounter;
    else
        OCMATAE.jumpcostatecoord=[];
        OCMATAE.jumpcostateindex=[];
    end
    if isempty(jump)
        jump=zeros(1,sccounter);
    end
    if length(jump)==1 && sccounter>1
        jump=repmat(jump,1,sccounter);
    end
elseif OCMATAE.stateconstraint && isempty(jump)
        OCMATAE.jumpcostatecoord=[];
        OCMATAE.jumpcostateindex=[];
end
numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1) ];
elseif length(sol.parameters)==numswitchtimes
    sol.parameters=[jump(:).' sol.parameters];
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1)];
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(1).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(1).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(1).odedim+(1:domaindata(1).aedim);
end
numode=domaindata(arcarg2arcindex(arcargument(ocAsym))).odedim;

% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
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
depvar=dependentvar(ocAsym);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.excludecoordinate4ep=excludecoordinate4ep;

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=sccounter+[1:numswitchtimes];
OCMATAE.inftimetransformation=inftimetransformation(ocAsym);
if isempty(OCMATAE.inftimetransformation)
    OCMATAE.inftimetransformation=0;
end
switch pathtpe
    case {'s','sc','cs','sts','ws'}
        OCMATAE.truncationtime=sol.arcinterval(end);
    case {'u','uc','cu','wu','stu'}
        OCMATAE.truncationtime=-abs(sol.arcinterval(end));
        sol.arcinterval=-abs(sol.arcinterval(end));
end
OCMATAE.saddlepoint=dependentvar(limSet);
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
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=objectivevaluecalc;
if isempty(OCMATAE.switchtimecoord)
    OCMATAE.movinghorizoncoord=1;
end
OCMATAE.movinghorizon=movinghorizon;
continuationparameter=0;
if ~isempty(freevector)
    OCMATAE.freevectorindex=length(sol.parameters)+(1:size(freevector,2));
    %sol.parameters=[sol.parameters ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex)];
    sol.parameters=[sol.parameters 0];
    OCMATAE.startvalue=OCMATAE.startvalue;%-ocAsym.solverinfo.freevector*ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex);
    OCMATAE.freevector=freevector(:);%ocAsym.solverinfo.freevector;
    continuationparameter=0;
    %OCMATAE.continuationvector=ocAsym.solverinfo.continuationvector;
else
    OCMATAE.freevector=[];
end
if movinghorizon
    yend=depvar(:,end);
    hatx=OCMATAE.saddlepoint;
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
        hatx(OCMATAE.excludecoordinate4ep,:)=[];
    end
    OCMATAE.distance=norm(yend(1:length(hatx))-hatx);
    OCMATAE.movinghorizoncoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
sol.parameters=[sol.parameters continuationparameter];
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocAsym,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocAsym,1)))];
end
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=size(sol.y,1);
end
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;

if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATAE.exogenousnumberofstates;
else
    OCMATAE.exogenousdynamicscoordinate=[];
end

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATAE.ODEcoord=1:numberofodes+OCMATAE.exogenousnumberofstates;

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
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