function sol=initocmat_AE_FTE_IS(mmObj,ocTrjMP,varargin)
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
freeparameter=0;
parindex=[];
targetcoordinate=[];
objectivevaluecalc=[];
userbc=[];
% input argument ocTrjMP is either a cell of octrajectories or a multi path object

ocTrjMP=ocmultipath(ocTrjMP);
indifforder=multiplicity(ocTrjMP);
nummod=numberofmodels(mmObj);

OCMATINDIF.infinitetime=zeros(1,indifforder);
for ii=1:indifforder
    ocTrjStruct=struct(ocTrjMP(ii));
    if length(ocTrjMP(ii).y(:,1))==2*statenum(mmObj.Model{1})
        
        o=objectivefunction(mmObj.Model{ii},ocTrjMP(ii),1);
        ocTrjStruct.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(mmObj.Model{ii},ocTrjMP(ii),1)))];
        ocTrjTmp{ii}=octrajectory(ocTrjStruct);
        
    else
        ocTrjTmp{ii}=ocTrjMP(ii);
    end
end
ocTrjMP=ocmultipath(ocTrjTmp);

if isempty(mmObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrjMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
userbcidx=find(strcmpi(varargin,'userbc'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
targettypeidx=find(strcmpi(varargin,'targettype'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
freevectoridx=find(strcmpi(varargin,'freevector'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));

if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end

if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
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
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if isempty(targettype)
    targettype='initialstate';
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=1;
end

scoord=statecoord(mmObj);
cscoord=costatecoord(mmObj);
for ii=1:nummod
    OCMATINDIF.statecostatecoord{ii}=[scoord{ii}(:).' cscoord{ii}(:).'];
end
switch lower(targettype)
    case 'initialstate'
        OCMATINDIF.targettype=1;
    case 'endtime'
        OCMATINDIF.targettype=2;
    case 'parameter'
        OCMATINDIF.targettype=3;
        parindex=varargin{targettypeidx+2};
        if ischar(parindex)
            parindex=parameterindex(mmObj,parindex);
        end
    otherwise
        OCMATINDIF.targettype=1;
end
if freeparameter
    parindex=varargin{freeparameteridx+1};
    if ischar(parindex)
        parindex=parameterindex(mmObj,parindex);
    end
    OCMATINDIF.freeparameter=1;
else
    OCMATINDIF.freeparameter=0;
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.objectivevaluecalc=1;
OCMATINDIF.hitvalue=hitvalue;
OCMATINDIF.userbcvalue=userbc;
% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     ocTrjMP(ii).solver='';
%     ocTrjMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=submodelname(mmObj);
OCMATCONT.modelfunc=modelspecificfunc(mmObj,'4IndifferenceSolutionContinuation');
OCMATINDIF.numberofmodels=nummod;
% initialize global variable (OCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
for ii=1:nummod
    funch{ii}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
end
for ii=1:nummod
    OCMATINDIF.canonicalsystem{ii}=funch{ii}{1};
    OCMATINDIF.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATINDIF.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};
    OCMATINDIF.canonicalsystemhessian{ii}=funch{ii}{3}{1};
    OCMATINDIF.canonicalsystemparameterhessian{ii}=funch{ii}{3}{2};

    % function for the boundary conditions
    OCMATINDIF.bcinitial{ii}=funch{ii}{5}{1};
    OCMATINDIF.bcasymptotic{ii}=funch{ii}{5}{2};
    OCMATINDIF.bctransversality{ii}=funch{ii}{5}{3};
    OCMATINDIF.bcindifference{ii}=funch{ii}{5}{5};
    OCMATINDIF.salvagevalue{ii}=funch{ii}{5}{6};
    OCMATINDIF.objectivevalue{ii}=funch{ii}{5}{7};
    if ~isempty(hitvalue)
        OCMATINDIF.hitvaluefunc{ii}=funch{ii}{5}{4};
    end
    if ~isempty(userbc)
        OCMATINDIF.userbc{ii}=funch{ii}{5}{8};
    else
        OCMATINDIF.userbc=[];
    end
    % function for Jacobian
    OCMATINDIF.bcjacobianinitial{ii}=funch{ii}{6}{1};
    OCMATINDIF.bcjacobiantransversality{ii}=funch{ii}{6}{2};

    % function describing the hybrid structure of the problem
    OCMATINDIF.hybridinfo{ii}=funch{ii}{7}{1};
    OCMATINDIF.domain{ii}=funch{ii}{7}{2};
    OCMATINDIF.guard{ii}=funch{ii}{7}{3};
    OCMATINDIF.reset{ii}=funch{ii}{7}{4};
    OCMATINDIF.switchtime{ii}=funch{ii}{7}{5};
    OCMATINDIF.jacobianguard{ii}=funch{ii}{7}{7};
    OCMATINDIF.jacobianreset{ii}=funch{ii}{7}{8};
    OCMATINDIF.domaindiscretization{ii}=funch{ii}{7}{9};

    OCMATINDIF.objectivefunction{ii}=funch{ii}{8}{1};
    OCMATINDIF.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
    OCMATINDIF.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
    OCMATINDIF.objectivefunctionderivativetime{ii}=funch{ii}{8}{4};

    if ~isautonomous(mmObj.Model{ii})
        OCMATINDIF.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
    end
    % general function
    OCMATINDIF.plotcontinuation{ii}=funch{ii}{11};
    OCMATINDIF.testadmissibility{ii}=funch{ii}{12};
    OCMATINDIF.datapath{ii}=funch{ii}{20};
    OCMATINDIF.saveintermediatefiles{ii}=funch{ii}{21};

    hybridinfo{ii}=OCMATINDIF.hybridinfo{ii}();
    
    for jj=1:numel(hybridinfo{ii}.arcarg)
        domaindata{ii}(jj)=OCMATINDIF.domain{ii}(hybridinfo{ii}.arcarg(jj));
    end
    for jj=1:numel(domaindata{ii})
        OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).numae=domaindata{ii}(jj).aedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).daeorder=domaindata{ii}(jj).daeorder;
        OCMATCONT.DOMAINDDATA{ii}(jj).numeq=numel(domaindata{ii}(jj).daeorder);%number of equations
        if objectivevaluecalc
            OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(1).odedim+1;
            if jj==1
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq+1;%number of equations
            else
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq;%number of equations
            end
        end
        OCMATCONT.DOMAINDDATA{ii}(jj).eqcoord=1:OCMATCONT.DOMAINDDATA{ii}(jj).numeq;
        OCMATCONT.DOMAINDDATA{ii}(jj).odecoord=1:domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).aecoord=domaindata{ii}(jj).odedim+(1:domaindata{ii}(ii).aedim);
    end
end

OCMATINDIF.parametervalue=parametervalue(mmObj);
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
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue{1}(parindex{1})];
    OCMATINDIF.parametercoord=length(sol.parameters)-length(parindex{1})+(1:length(parindex{1}));
    OCMATINDIF.endtime=arcintv(end);
end
if  OCMATINDIF.freeparameter
    %OCMATINDIF.freevectorcoord=[];
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue{1}(parindex{1})];
    OCMATINDIF.parametercoord=length(sol.parameters)-length(parindex{1})+(1:length(parindex{1}));
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
OCMATINDIF.startvalue=depvar(scoord{1},1);

pathname=OCMATINDIF.datapath{1}();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles{1}();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.autonomous=isautonomous(mmObj);
if objectivevaluecalc
    OCMATINDIF.objectivevaluecoord=size(sol.y,1);
end
if ~objectivevaluecalc && length(sol.y(:,1))>2*statenum(mmObj.Model{1})
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
