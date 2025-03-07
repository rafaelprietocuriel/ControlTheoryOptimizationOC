function sol=initocmat_FIMP_T_IS(ocObj,hocMP,varargin)
%
% INITOCMAT_AE_IS initialization for the continuation of an indifference
% threshold
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE) the
% continuation of an indifference threshold in the state space is
% initialized.
% OCOBJ          ... corresponding optimal control model
% HOCMP           ... two cell array of hybridoctrajectories, for the different
%                    solution paths or an instance of an hybridocmultipath
%                    object.
% TARGETCOORDINATE ... the continuation is done along time horizon
% TARGETVALUE    ... The value of the target vector for the time horizon.
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and IOCMATINDIF
% are initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SOL=INITOCMAT_FIMP_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

clear global OCMATCONT IOCMATINDIF

global OCMATCONT IOCMATINDIF
sol=[];
targetvalue=[];
targettype='';
freevector=[];

% input argument hocMP is either a cell of ocasymptotics or a multi path object
hocMP=hybridocmultipath(hocMP);
indifforder=multiplicity(hocMP);
for ii=1:indifforder
    pathtpe{ii}=pathtype(hocMP(ii));
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(hocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targettypeidx=find(strcmpi(varargin,'targettype'));
freevectoridx=find(strcmpi(varargin,'freevector'));

if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end

snum=statenum(ocObj);
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
IOCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

switch lower(targettype)
    case 'initialstate'
        IOCMATINDIF.targettype=1;
    case 'endtime'
        IOCMATINDIF.targettype=2;
    case 'parameter'
        IOCMATINDIF.targettype=3;
    otherwise
        IOCMATINDIF.targettype=1;
end
IOCMATINDIF.freevector=freevector;
IOCMATINDIF.objectivevaluecalc=1;
% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     hocMP(ii).solver='';
%     hocMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4IndifferenceSolutionContinuation');

% initialize global variable (IOCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for indifference solution continuation

IOCMATINDIF.canonicalsystem=funch{1};
IOCMATINDIF.canonicalsystemjacobian=funch{2}{1};
IOCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
% IOCMATINDIF.canonicalsystemhessian=funch{3}{1};
% IOCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
IOCMATINDIF.bcinitial=funch{5}{1};
IOCMATINDIF.bctransversality=funch{5}{2};
IOCMATINDIF.bcevent=funch{5}{3};
IOCMATINDIF.bcinteriorevent=funch{5}{4};
IOCMATINDIF.bcindifference=funch{5}{5};

% function for Jacobian
IOCMATINDIF.bcjacobianinitial=funch{6}{1};
IOCMATINDIF.bcoptimalhorizon=funch{5}{3};
IOCMATINDIF.bcevent=funch{5}{4};
IOCMATINDIF.bcinteriorevent=funch{5}{5};

% function describing the hybrid structure of the problem
IOCMATINDIF.hybridinfo=funch{7}{1};
IOCMATINDIF.domain=funch{7}{2};
IOCMATINDIF.guard=funch{7}{3};
IOCMATINDIF.reset=funch{7}{4};
IOCMATINDIF.switchtime=funch{7}{5};
IOCMATINDIF.jacobianguard=funch{7}{7};
IOCMATINDIF.jacobianreset=funch{7}{8};
IOCMATINDIF.domaindiscretization=funch{7}{9};
IOCMATINDIF.objectivefunction=funch{8}{1};
IOCMATINDIF.objectivefunctionjacobian=funch{8}{2};
IOCMATINDIF.objectivefunctionparameterjacobian=funch{8}{3};
if ~isautonomous(ocObj)
    IOCMATINDIF.canonicalsystemderivativetime=funch{2}{3};
    IOCMATINDIF.objectivefunctionderivativetime=funch{8}{4};
end
IOCMATINDIF.objectivevalue=funch{8}{5};

% general function
IOCMATINDIF.plotcontinuation=funch{11};
IOCMATINDIF.testadmissibility=funch{12};
IOCMATINDIF.datapath=funch{20};
IOCMATINDIF.saveintermediatefiles=funch{21};


hybridinfo=IOCMATINDIF.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=IOCMATINDIF.domain(hybridinfo.arcarg(ii));
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(1).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(1).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim+1;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder)+1;%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(1).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(1).odedim+(1:domaindata(ii).aedim);
end

sol=generatesolstruct(hocMP);

% mode and path specific variables
OCMATCONT.HE.edge=[];
IOCMATINDIF.jumparg=[];
arcoffset=0;
paroffset=0;
jumpoffset=0;
for ii=1:indifforder
    arcn=arcnum(hocMP(ii));
    depvar=dependentvar(hocMP(ii));

    arctime=arcinterval(hocMP(ii));
    IOCMATINDIF.arctime{ii}=arctime;
    jumparg=jumpargument(hocMP(ii));
    IOCMATINDIF.jumparg=[IOCMATINDIF.jumparg jumparg];

    switchtimeidx=find(jumparg>0 | [-1 jumparg(2:end-1) -1]==0);
    IOCMATINDIF.initialdepvarcoord{ii}=[scoord(:).' cscoord(:).']+paroffset;
    IOCMATINDIF.enddepvarcoord{ii}=2*snum+[scoord(:).' cscoord(:).']+paroffset;
    IOCMATINDIF.switchtimecoord{ii}=4*snum+(1:length(switchtimeidx))+paroffset;
    IOCMATINDIF.switchtimeidx{ii}=switchtimeidx;
    IOCMATINDIF.varyarcintervalidx{ii}=zeros(1,length(arctime));
    IOCMATINDIF.varyarcintervalidx{ii}(switchtimeidx)=IOCMATINDIF.switchtimecoord{ii};
    arcarg=arcargument(hocMP(ii));
    IOCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    IOCMATINDIF.numarc(ii)=arcn;
    IOCMATINDIF.initialstateindex(ii)=numel(hocMP(ii).x);
    IOCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    IOCMATINDIF.jumpcoord{ii}=(1:arcn+1)+jumpoffset;
    arcoffset=arcoffset+arcn;
    jumpoffset=jumpoffset+arcn+1;
    sol.parameters=[sol.parameters depvar([scoord cscoord],1).' depvar([scoord cscoord],end).'];
    sol.parameters=[sol.parameters arctime(switchtimeidx)];
    paroffset=paroffset+length(sol.parameters);
end

if IOCMATINDIF.targettype==1 % initial state continuation
    sol.parameters=[sol.parameters arctime(end)];
    IOCMATINDIF.endtimecoord=length(sol.parameters);
    IOCMATINDIF.freevectorcoord=length(sol.parameters)+1;
    IOCMATINDIF.targetvalue=targetvalue;
    IOCMATINDIF.freevector=targetvalue-depvar(scoord,1);
    sol.parameters=[sol.parameters 0];
elseif IOCMATINDIF.targettype==2 % end time continuation
    if ~isempty(freevector)
        sol.parameters=[sol.parameters 0]; % more than one states
        IOCMATINDIF.freevectorcoord=length(sol.parameters);
    else
        IOCMATINDIF.freevectorcoord=[];
    end
    sol.parameters=[sol.parameters  arctime(end)];
    IOCMATINDIF.endtimecoord=length(sol.parameters);
    IOCMATINDIF.targetvalue=targetvalue;
else
end
for ii=1:indifforder
    IOCMATINDIF.varyarcintervalidx{ii}(end)=IOCMATINDIF.endtimecoord;
end

IOCMATINDIF.initialstateindex=cumsum(IOCMATINDIF.initialstateindex);
IOCMATINDIF.initialstateindex=[0 IOCMATINDIF.initialstateindex(end-1)]+1;
IOCMATINDIF.solutionindex=zeros(1,sum(IOCMATINDIF.numarc));% relate arcindex to indifference solution
counter=0;
for order=1:indifforder
    counter_start=counter+1;
    counter=counter+IOCMATINDIF.numarc(order);
    IOCMATINDIF.solutionindex(counter_start:counter)=order;
end

IOCMATINDIF.cumsumnumarc=cumsum(IOCMATINDIF.numarc);
IOCMATINDIF.initcoord=[1 IOCMATINDIF.cumsumnumarc(1:end-1)+1];

IOCMATINDIF.indifferenceorder=indifforder;
IOCMATINDIF.parametervalue=parametervalue(ocObj);
IOCMATINDIF.initialtime=sol.x0;

IOCMATINDIF.statecoordinate=scoord;
IOCMATINDIF.initialstate=depvar(scoord,1);

IOCMATINDIF.pathtype=pathtpe;
pathname=IOCMATINDIF.datapath();
[resultfile,globalvarfile]=IOCMATINDIF.saveintermediatefiles();
IOCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
IOCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
IOCMATINDIF.autonomous=isautonomous(ocObj);
IOCMATINDIF.objectivevaluecoord=size(sol.y,1);

OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;
OCMATCONT.continuation=1;

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.x([1 end])=[];
sol.y=dependentvar(ocMultiPath(1));
sol.y(:,[1 end])=[];
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    x=independentvar(ocMultiPath(ii))-x0+sol.x(end);
    x([1 end])=[];
    sol.x=[sol.x x];
    y=dependentvar(ocMultiPath(ii));
    y(:,[1 end])=[];
    sol.y=[sol.y y];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
end
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
sol.parameters=[];
