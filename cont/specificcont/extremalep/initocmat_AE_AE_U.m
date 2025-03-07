function sol=initocmat_AE_AE_U(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
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
modelfuncidx=find(strcmpi(varargin,'modelfunction'));
if ~isempty(jumpidx)
    jump=varargin{jumpidx+1};
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if isempty(movinghorizon)
    movinghorizon=0;
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if isempty(modelfuncidx)
    modelfunc='4SaddlePathContinuation';
else
    modelfunc=varargin{modelfuncidx+1};
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,modelfunc);

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

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianasymptotic=funch{6}{2};
OCMATAE.bcjacobiantransversality=funch{6}{3};


OCMATAE.objectivefunction=funch{8}{1};
OCMATAE.objectivefunctionjacobian=funch{8}{2};
OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};

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

try
    OCMATAE.solutionadapt=funch{50};
catch
    OCMATAE.solutionadapt=[];
end
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
J=linearization(limSet);
depvar=dependentvar(ocAsym);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=sccounter+[1:numswitchtimes];
OCMATAE.inftimetransformation=inftimetransformation(ocAsym);
if isempty(OCMATAE.inftimetransformation)
    OCMATAE.inftimetransformation=0;
end
switch pathtpe
    case {'s','sc','cs'}
        OCMATAE.truncationtime=sol.arcinterval(end);
    case 'u'
        OCMATAE.truncationtime=-sol.arcinterval(end);
end
OCMATAE.linearization=J;
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=false;
if isempty(OCMATAE.switchtimecoord)
    OCMATAE.movinghorizoncoord=1;
end
OCMATAE.movinghorizon=movinghorizon;
continuationparameter=0;
if ~isempty(freevector)
    OCMATAE.freevectorindex=length(sol.parameters)+(1:size(freevector,2));
    sol.parameters=[sol.parameters ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex)];
    OCMATAE.startvalue=OCMATAE.startvalue-ocAsym.solverinfo.freevector*ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex);
    OCMATAE.freevector=ocAsym.solverinfo.freevector;
    continuationparameter=0;
    OCMATAE.continuationvector=ocAsym.solverinfo.continuationvector;
else
    OCMATAE.freevector=[];
end
if movinghorizon
    OCMATAE.distance=norm(OCMATAE.saddlepoint-depvar(:,end));
    OCMATAE.movinghorizoncoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
sol.parameters=[sol.parameters continuationparameter];

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;

if ~isempty(OCMATAE.solutionadapt)
    sol=OCMATAE.solutionadapt(sol);
end

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