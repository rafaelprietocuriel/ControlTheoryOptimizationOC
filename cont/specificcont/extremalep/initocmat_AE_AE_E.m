function sol=initocmat_AE_AE_E(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
%
% INITOCMAT_AE_AE_E initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_AE_E(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
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
endpoint=[]; % provide possible jumps of costate for pure state constraints
pathtpe=pathtype(ocAsym);
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
if nargin>=6
    endpoint=varargin{1};
end
if nargin>=7
    jump=varargin{2};
end
targetvalue=targetvalue(:);

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
OCMATAE.bcexactasymptotic=funch{5}{2};

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianexactasymptotic=funch{6}{2};

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

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
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
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
depvar=dependentvar(ocAsym);
if isempty(endpoint)
    endpoint=depvar(:,end);
end

OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=initialtime(ocAsym);
OCMATAE.endarcarg=arcargument(limSet);

OCMATAE.inftimetransformation=0;
OCMATAE.linearization=J;
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;

sol=generatesolstruct(ocAsym,ocObj,endpoint,jump);
sol.parameters=[sol.parameters 0];

function sol=generatesolstruct(ocTrj,ocObj,endpoint,jump)
global OCMATAE

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=parameters(ocTrj);
solverinfo0=solverinfo(ocTrj);
arcarg=[sol.arcarg OCMATAE.endarcarg];

% test for state constraints
sccounter=0;
if OCMATAE.stateconstraint
    arcarg4sc=arcargwithactivestateconstraint(ocObj);
    for ii=2:length(arcarg)
        if ismember(arcarg(ii),arcarg4sc)
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
end
sol.parameters=[endpoint(:).' sol.arcinterval(2:end)];

OCMATAE.exactendpointcoord=1:length(endpoint);
OCMATAE.switchtimecoord=OCMATAE.exactendpointcoord(end)+[1:length(sol.arcinterval(2:end))];
if ~isempty(jump)
    OCMATAE.jumpcostatecoord=solverinfo0.switchtimecoord(end)+[1:sccounter];
    sol.parameters=[sol.parameters jump];
end
sol.solver='';

if isfield(solverinfo0,'tangent')
    sol.solverinfo.tangent=solverinfo0.tangent;
else
    sol.solverinfo.tangent=[];
end
if isfield(solverinfo0,'coeff')
    sol.solverinfo.coeff=solverinfo0.coeff;
else
    sol.solverinfo.coeff=[];
end
if isfield(solverinfo0,'tmesh')
    sol.solverinfo.tmesh=solverinfo0.tmesh;
else
    sol.solverinfo.coeff=[];
end
if isfield(solverinfo0,'yp')
    sol.solverinfo.yp=solverinfo0.yp;
end
if isfield(solverinfo0,'ypmid')
    sol.solverinfo.ypmid=solverinfo0.ypmid;
end
if isfield(solverinfo0,'ymid')
    sol.solverinfo.yp=solverinfo0.yp;
end
if isfield(solverinfo0,'ymid')
    sol.solverinfo.ymid=solverinfo0.ymid;
end

OCMATAE.arcarg=arcarg;
OCMATAE.numarc=length(arcarg);
OCMATAE.edge=[arcarg(1:end-1); ...
    arcarg(2:end)];
