function sol=initocmat_AE_FTE_IS(ocObj,ocTrjMP,varargin)
%
% initocmat_FTE_IS initialization for the continuation of an indifference
% threshold for a mixed finite and infinite time horizon problem

clear global OCMATCONT OCMATINDIF

global OCMATCONT OCMATINDIF
sol=[];
targetvalue=[];
hitvalue=[];
targettype='';
freevector=[];
parindex=[];
targetcoordinate=[];
objectivevaluecalc=[];
% input argument ocTrjMP is either a cell of octrajectories or a multi path object

ocTrjMP=ocmultipath(ocTrjMP);
indifforder=multiplicity(ocTrjMP);
OCMATINDIF.infinitetime=zeros(1,indifforder);
for ii=1:indifforder
    if isocasymptotic(ocTrjMP(ii))
        OCMATINDIF.infinitetime(ii)=1;
    end
    ocTrjStruct=struct(ocTrjMP(ii));
    if ~isocasymptotic(ocTrjMP(ii)) && length(ocTrjMP(ii).y(:,1))==2*statenum(ocObj)
        
        o=objectivefunction(ocObj,ocTrjMP(ii),1);
        ocTrjStruct.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocTrjMP(ii),1)))];
        ocTrjTmp{ii}=octrajectory(ocTrjStruct);
        
    else
        ocTrjTmp{ii}=ocTrjMP(ii);
    end
end
ocTrjMP=ocmultipath(ocTrjTmp);

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrjMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
targettypeidx=find(strcmpi(varargin,'targettype'));
freevectoridx=find(strcmpi(varargin,'freevector'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));

if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end

if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if isempty(targettype)
    targettype='initialstate';
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end

scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

switch lower(targettype)
    case 'initialstate'
        OCMATINDIF.targettype=1;
    case 'endtime'
        OCMATINDIF.targettype=2;
    case 'parameter'
        OCMATINDIF.targettype=3;
        parindex=varargin{targettypeidx+2};
        if ischar(parindex)
            parindex=parameterindex(ocObj,parindex);
        end
    otherwise
        OCMATINDIF.targettype=1;
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.objectivevaluecalc=1;
OCMATINDIF.hitvalue=hitvalue;

% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     ocTrjMP(ii).solver='';
%     ocTrjMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4IndifferenceSolutionContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for indifference solution continuation

OCMATINDIF.canonicalsystem=funch{1};
OCMATINDIF.canonicalsystemjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
% OCMATINDIF.canonicalsystemhessian=funch{3}{1};
% OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.bcindifference=funch{5}{5};
OCMATINDIF.salvagevalue=funch{5}{6};
OCMATINDIF.objectivevalue=funch{5}{7};

% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcoptimalhorizon=funch{5}{3};
OCMATINDIF.bcevent=funch{5}{4};
OCMATINDIF.bcinteriorevent=funch{5}{5};

% function describing the hybrid structure of the problem
OCMATINDIF.hybridinfo=funch{7}{1};
OCMATINDIF.domain=funch{7}{2};
OCMATINDIF.guard=funch{7}{3};
OCMATINDIF.reset=funch{7}{4};
OCMATINDIF.switchtime=funch{7}{5};
OCMATINDIF.jacobianguard=funch{7}{7};
OCMATINDIF.jacobianreset=funch{7}{8};
OCMATINDIF.domaindiscretization=funch{7}{9};
OCMATINDIF.objectivefunction=funch{8}{1};
OCMATINDIF.objectivefunctionjacobian=funch{8}{2};
if objectivevaluecalc
    OCMATINDIF.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATINDIF.objectivefunctionderivativetime=funch{8}{4};
end

if ~isautonomous(ocObj)
    OCMATINDIF.canonicalsystemderivativetime=funch{2}{3};
end

% general function
OCMATINDIF.plotcontinuation=funch{11};
OCMATINDIF.testadmissibility=funch{12};
OCMATINDIF.datapath=funch{20};
OCMATINDIF.saveintermediatefiles=funch{21};


hybridinfo=OCMATINDIF.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATINDIF.domain(hybridinfo.arcarg(ii));
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(1).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(1).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(1).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(1).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(1).odedim+(1:domaindata(ii).aedim);
end
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.objectivevaluecalc=objectivevaluecalc;

sol=generatesolstruct(ocTrjMP);

OCMATINDIF.parameterindex=parindex;
% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
for ii=1:indifforder
    arcn=arcnum(ocTrjMP(ii));
    arcintv=arcinterval(ocTrjMP(ii));
    arcarg=arcargument(ocTrjMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocTrjMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcoffset+arcn;
    sol.parameters=[sol.parameters switchtimes{ii}];
end
counter=0;
paroffset=0;
for ii=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    if ii==1
        if ~isempty(switchtimes{ii})
            OCMATINDIF.switchtimecoord{ii}=(1:length(switchtimes{ii}));
        else
            OCMATINDIF.switchtimecoord{ii}=[];
        end
    else
        OCMATINDIF.switchtimecoord{ii}=paroffset+(1:length(switchtimes{ii}));
    end
    if ~isempty(OCMATINDIF.switchtimecoord{ii})
        paroffset=OCMATINDIF.switchtimecoord{ii}(end);
    end
end
depvar=dependentvar(ocTrjMP(1));

if OCMATINDIF.targettype==1 % initial state continuation
    OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
    for ii=1:size(freevector,2)
        sol.parameters=[sol.parameters 0];
    end
    OCMATINDIF.freevector=freevector;
    OCMATINDIF.endtime=arcintv(end);
elseif OCMATINDIF.targettype==2 % end time continuation
    if ~isempty(freevector)
        sol.parameters=[sol.parameters 0]; % more than one states
        OCMATINDIF.freevectorcoord=length(sol.parameters);
    else
        OCMATINDIF.freevectorcoord=[];
    end
    sol.parameters=[sol.parameters  arcintv(end)];
    OCMATINDIF.endtimecoord=length(sol.parameters);
elseif OCMATINDIF.targettype==3
    if ~isempty(freevector)
        sol.parameters=[sol.parameters 0]; % more than one states
        OCMATINDIF.freevectorcoord=length(sol.parameters);
    else
        OCMATINDIF.freevectorcoord=[];
    end
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];
    OCMATINDIF.parametercoord=length(sol.parameters)-length(parindex)+(1:length(parindex));
    OCMATINDIF.endtime=arcintv(end);
end
OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;

OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution
counter=0;
for order=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(order);
    OCMATINDIF.solutionindex(counter_start:counter)=order;
end

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.endcoord=OCMATINDIF.cumsumnumarc(1:end);

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.initialtime=sol.x0;

OCMATINDIF.statecoordinate=scoord;
OCMATINDIF.startvalue=depvar(scoord,1);

pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.autonomous=isautonomous(ocObj);
if objectivevaluecalc
    OCMATINDIF.objectivevaluecoord=size(sol.y,1);
end
if ~objectivevaluecalc && length(sol.y(:,1))>2*statenum(ocObj)
    sol.y(end,:)=[];
end

OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;
OCMATCONT.continuation=1;

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
sol.parameters=parameters(ocMultiPath(1));
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
end
sol.parameters=[];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
