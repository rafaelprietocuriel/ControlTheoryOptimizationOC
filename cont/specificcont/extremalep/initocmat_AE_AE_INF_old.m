function sol=initocmat_AE_AE_INF(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
%
% initocmat_AE_AE_INF initialization for asymptotic extremal calculation
%
% SOL=initocmat_AE_AE_INF(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
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
convergingcoord=[];
pathtpe=pathtype(ocAsym);
% transversalityvalue=[];

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
freehorizonidx=find(strcmpi(varargin,'freehorizon'),1);
% freevectoridx=find(strcmpi(varargin,'freevector'));
initialstateidx=find(strcmpi(varargin,'initialstate'),1);
% transversalityvalueidx=find(strcmpi(varargin,'transversalityvalue'));
convergingcoordidx=find(strcmpi(varargin,'convergingcoord'));
% if ~isempty(freevectoridx)
%     freevector=varargin{freevectoridx+1};
% end
% if ~isempty(transversalityvalueidx)
%     transversalityvalue=varargin{transversalityvalueidx+1};
% end
if ~isempty(convergingcoordidx)
    convergingcoord=varargin{convergingcoordidx+1};
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

if ~isempty(initialstateidx)
   conttype='initialstate';
elseif ~isempty(freehorizonidx)
    conttype='freehorizon';
% elseif ~isempty(transversalityvalueidx)
%     conttype='transversalityvalue';
% elseif ~isempty(freevectoridx)
%     conttype='freevector';
else
    ocmaterror('Type of continuation cannot be determined.')
end
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
OCMATAE.bcinftransversality=funch{5}{3};

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

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
sccounter=0;
numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=sol.arcinterval(2:end-1);
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
switch pathtpe
    case {'s','sc','cs'}
        OCMATAE.truncationtime=sol.arcinterval(end);
    case 'u'
        OCMATAE.truncationtime=-sol.arcinterval(end);
end
OCMATAE.linearization=J;
OCMATAE.saddlepoint=dependentvar(limSet);
if ~isempty(convergingcoord)
    Jred=J(convergingcoord,:);
    Jred=Jred(:,convergingcoord);
    OCMATAE.asymptoticmatrix=asymptoticbc(Jred,pathtpe,'c',ZeroDeviationTolerance);
else
    OCMATAE.asymptoticmatrix=[];
end
OCMATAE.convergingcoord=convergingcoord;
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.freevector=[];
OCMATAE.initialstate=depvar(statecoord(ocObj),1);
OCMATAE.statecoordinate=statecoord(ocObj);
OCMATAE.freehorizon=false;


switch conttype
    case 'initialstate'
        OCMATAE.startvalue=depvar(contcoordinate,1);
        OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
        OCMATAE.statecoordinate=contcoordinate;
        OCMATAE.initialstatecont=true;
        continuationparameter=0;
    case 'freehorizon'
        OCMATAE.startvalue=sol.arcinterval(contcoordinate);
        OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
        OCMATAE.freehorizon=true;
        continuationparameter=0;
        OCMATAE.initialstatecont=false;
%     case 'transversalityvalue'
%         OCMATAE.transversalitycoord=contcoordinate;
%         OCMATAE.startvalue=transversalityvalue;
%         OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
%         continuationparameter=0;
%     case 'freevector'
%         OCMATAE.freevectorindex=length(sol.parameters)+(1:size(freevector,2));
%         sol.parameters=[sol.parameters ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex)];
%         OCMATAE.startvalue=OCMATAE.startvalue-ocAsym.solverinfo.freevector*ocAsym.solverinfo.parameters(ocAsym.solverinfo.freevectorindex);
%         OCMATAE.freevector=ocAsym.solverinfo.freevector;
%         continuationparameter=0;
%         OCMATAE.continuationvector=ocAsym.solverinfo.continuationvector;
end
OCMATAE.truncationtime=sol.arcinterval(end);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.continuationindex=length(sol.parameters)+1;
OCMATAE.continuationspecification=conttype;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=false;
sol.parameters=[sol.parameters continuationparameter];

switch conttype
    case 'initialstate'
        OCMATAE.continuationcoord=length(sol.parameters);
        OCMATAE.freehorizoncoord=[];
    case 'freehorizon'
        OCMATAE.freehorizoncoord=length(sol.parameters);
end

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