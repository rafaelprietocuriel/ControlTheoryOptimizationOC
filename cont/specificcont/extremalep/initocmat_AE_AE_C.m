function sol=initocmat_AE_AE_C(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
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
if ~isempty(jumpidx)
    jump=varargin{jumpidx+1};
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if isempty(movinghorizon)
    movinghorizon=[];
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
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.bcinitialpointcurve=funch{5}{4};

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

depvar=dependentvar(ocAsym);
sol=generatesolstruct(ocAsym,solver(ocAsym));

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
sccounter=0;
if OCMATAE.stateconstraint
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
end
numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1) ];
elseif length(sol.parameters)==numswitchtimes
    sol.parameters=[jump(:).' sol.parameters];
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=[jump(:).' sol.arcinterval(2:end-1)];
end
sol.parameters=[sol.parameters depvar(contcoordinate)];
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
OCMATAE.targetvalue=targetvalue;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=false;
if isempty(OCMATAE.switchtimecoord)
    movinghorizon=[];
end
OCMATAE.movinghorizon=movinghorizon;

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;

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
        case 'extremalc2ep'
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