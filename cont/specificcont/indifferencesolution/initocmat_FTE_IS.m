function sol=initocmat_FTE_IS(ocObj,ocTrjMP,continuationtype,continuationindex,continuationtarget,varargin)
%
% initocmat_FTE_IS initialization for the continuation of an indifference
% threshold for a finite time horizon problem

clear global OCMATCONT OCMATINDIF

global OCMATCONT OCMATINDIF
sol=[];
targetvalue=[];
freestatevector=[];
freeparametervector=[];
freeparametertarget=[];
freetimevector=[];
fixinitstate=[];
parindex=[];
objectivevaluecalc=[];
exogenousfunction=[];
jumpargument=[];
entryindex=[];
hitvalue=[];
variationalcalculation=[];
vjumpargument=[];
vfreetime=[];
static=0; % y coordinate, where value for the static optimization is computed
staticparameterindex=[];
userbc=[];
% input argument ocTrjMP is either a cell of octrajectories or a multi path object

ocTrjMP=ocmultipath(ocTrjMP);
indifforder=multiplicity(ocTrjMP);
for ii=1:indifforder
    ocTrjStruct=struct(ocTrjMP(ii));
    if length(ocTrjMP(ii).y(:,1))==2*statenum(ocObj)
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
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
freestatevectoridx=find(strcmpi(varargin,'freestatevector'));
freetimevectoridx=find(strcmpi(varargin,'freetimevector'));
freeparametervectoridx=find(strcmpi(varargin,'freeparametervector'));
freeparametertargetidx=find(strcmpi(varargin,'freeparametertarget'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
jumpididx=find(strcmpi(varargin,'jumpid'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
variationalcalculationidx=find(strcmpi(varargin,'variationalcalculation'));
vjumpargumentidx=find(strcmpi(varargin,'variationaljumpargument'));
vfreetimeidx=find(strcmpi(varargin,'variationalfreetime'));
staticidx=find(strcmpi(varargin,'static'));
staticparameteridx=find(strcmpi(varargin,'staticparameter'));
userbcidx=find(strcmpi(varargin,'userbc'));

if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(freestatevectoridx)
    freestatevector=varargin{freestatevectoridx+1};
end
if ~isempty(freeparametervectoridx)
    freeparametervector=varargin{freeparametervectoridx+1};
end
if ~isempty(freeparametertargetidx)
    freeparametertarget=varargin{freeparametertargetidx+1};
end
if ~isempty(freetimevectoridx)
    freetimevector=varargin{freetimevectoridx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end

if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
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
if ~isempty(jumpargumentidx)
    jumpargument=varargin{jumpargumentidx+1};
end
if ~isempty(entryindexidx)
    entryindex=varargin{entryindexidx+1};
end
if ~isempty(jumpididx)
    jumpid=varargin{jumpididx+1};
end
if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if ~isempty(variationalcalculationidx)
    variationalcalculation=varargin{variationalcalculationidx+1};
end
if ~isempty(vjumpargumentidx)
    vjumpargument=varargin{vjumpargumentidx+1};
end
if ~isempty(vfreetimeidx)
    vfreetime=varargin{vfreetimeidx+1};
end
if ~isempty(staticidx)
    static=varargin{staticidx+1};
end
if ~isempty(staticparameteridx)
    staticparameter=varargin{staticparameteridx+1};
    if ischar(staticparameter)
        staticparameterindex=parameterindex(ocObj,staticparameter);
    else
        staticparameterindex=staticparameter;
    end
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if isempty(static)
    static=0;
end

OCMATINDIF.variationalcalculation=variationalcalculation;
OCMATINDIF.static=static; % equals the coordinate for which the static optimization is carried out
OCMATINDIF.userbc=userbc;
for ii=1:indifforder
    if iscell(jumpargument)
        OCMATINDIF.stateconstraint{ii}=stateconstraint(ocObj)&&~isempty(jumpargument{ii});
    else
        OCMATINDIF.stateconstraint{ii}=0;
    end
end
if variationalcalculation
    OCMATINDIF.variationalinitialstates=zeros(2*statenum(ocObj),1);
    OCMATINDIF.variationalnumberofstates=length(OCMATINDIF.variationalinitialstates);
end
OCMATINDIF.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end

scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

switch lower(continuationtype)
    case 'initialstate'
        OCMATINDIF.continuationtype=0;
    case 'endtime'
        OCMATINDIF.continuationtype=1;
    case 'parameter'
        OCMATINDIF.continuationtype=2;
        if ischar(continuationindex)
            parindex=parameterindex(ocObj,continuationindex);
        end
end
OCMATINDIF.fixinitstatecoord=fixinitstate;
OCMATINDIF.freestatevector=freestatevector;
OCMATINDIF.freeparametervector=freeparametervector;
if ~isempty(freeparametervector)
    OCMATINDIF.freeparameterindex=parameterindex(ocObj,freeparametervector);
end
OCMATINDIF.freeparametertarget=freeparametertarget;
if ~isempty(freeparametertarget)
    OCMATINDIF.freeparametertarget=parameterindex(ocObj,freeparametervector);
    OCMATINDIF.freeparametertarget=OCMATINDIF.freeparametertarget(:).';
end
if ~isempty(freeparametertarget) && (length(OCMATINDIF.freeparameterindex)>1 && length(OCMATINDIF.freeparameterindex)~=length(freeparametertarget))
    % one dimesion of freedom is assumed
    ocmatmsg('Number of free parameter values and size of direction are not consistent')
    return
end
OCMATINDIF.freetimevector=freetimevector;

OCMATINDIF.initialstate=state(ocObj,ocTrjMP(1));
OCMATINDIF.initialstate=OCMATINDIF.initialstate(:,1);
OCMATINDIF.continuationtarget=continuationtarget(:);
OCMATINDIF.objectivevaluecalc=1;
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

if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamics=funch{4}{1};
    OCMATINDIF.exogenousjacobian=funch{4}{2};
    OCMATINDIF.exogenousparameterjacobian=funch{4}{3};
    OCMATINDIF.exogenousinitialstates=funch{4}{11};
end
if OCMATINDIF.variationalcalculation
    OCMATINDIF.variationaldynamics=funch{4}{4};
    OCMATINDIF.variationaljacobian=funch{4}{5};
    OCMATINDIF.variationalparameterjacobian=funch{4}{6};
    OCMATINDIF.variationalhamiltonian=funch{4}{10};
    OCMATINDIF.variationalguard=funch{5}{9};
    OCMATINDIF.variationalreset=funch{5}{10};
    OCMATINDIF.variationalbcinitial=funch{5}{11};
    OCMATINDIF.variationalbctransversality=funch{5}{12};
    if any([OCMATINDIF.stateconstraint{:}])
        OCMATINDIF.variationalbcstateconstraint=funch{5}{13};
        OCMATINDIF.variationalbctransversalitysc=funch{5}{14};
    end
    if OCMATINDIF.exogenousfunction
        OCMATINDIF.exogenousvariationaldynamics=funch{4}{7};
        OCMATINDIF.exogenousjacobian4variationalargument=funch{4}{8};
        OCMATINDIF.exogenousvariationaldynamicsjacobian=funch{4}{9};
        %OCMATINDIF.exogenousinitialstates=funch{4}{11};
    end
    OCMATINDIF.variationalobjectivefunction=funch{8}{5};
    OCMATINDIF.variationalobjectivefunctionjacobian=funch{8}{6};
    OCMATINDIF.variationalobjectivefunctionparameterjacobian=funch{8}{7};
    OCMATINDIF.variationalobjectivefunctionderivativetime=funch{8}{8};
    OCMATINDIF.variationalsalvagevalue=funch{8}{9};
end

% OCMATINDIF.canonicalsystemhessian=funch{3}{1};
% OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.bcindifference=funch{5}{5};
OCMATINDIF.salvagevalue=funch{5}{6};
if any([OCMATINDIF.stateconstraint{:}])
    OCMATINDIF.bcstateconstraint=funch{5}{7};
    OCMATINDIF.bctransversalitysc=funch{5}{8};
    OCMATINDIF.testjumpargument=funch{13};
end
if ~isempty(hitvalue)
    OCMATINDIF.hitvaluefunc=funch{5}{4};
end
if ~isempty(userbc)
    OCMATINDIF.userbcfunc=funch{5}{15};
end

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
OCMATINDIF.objectivefunctionparameterjacobian=funch{8}{3};
OCMATINDIF.objectivefunctionderivativetime=funch{8}{4};

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
OCMATINDIF.hitvalue=hitvalue;
%OCMATINDIF.modelparameter=modelparameter(ocTrjMP);
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
    freeparametercoordinatep1=length(sol.parameters)+1;
    sol.parameters=[sol.parameters switchtimes{ii}];
    OCMATINDIF.switchtimecoord{ii}=freeparametercoordinatep1:length(sol.parameters);
    if variationalcalculation
        %         if length(vfreetime)~=length(OCMATINDIF.freeswitchingtimeindex)
        %             ocmaterror('Number of variational free time arguments and number of free time arguments are different.')
        %         end
        freeparametercoordinatep1=length(sol.parameters)+1;
        sol.parameters=[sol.parameters vfreetime{ii}];
        OCMATINDIF.vfreetimecoord{ii}=freeparametercoordinatep1:length(sol.parameters);
    end
end
counter=0;
for ii=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
end
for ii=1:indifforder
    if OCMATINDIF.stateconstraint{ii}
        arcintv=arcinterval(ocTrjMP(ii));
        OCMATINDIF.jumpargument{ii}=jumpargument{ii};
        OCMATINDIF.entryindex{ii}=zeros(1,length(arcintv));
        OCMATINDIF.jumpid{ii}=zeros(1,length(arcintv));
        OCMATINDIF.entryindex{ii}(entryindex{ii})=1:length(jumpargument{ii});
        OCMATINDIF.jumpid{ii}(entryindex{ii})=jumpid{ii};
        if ~isempty(jumpargument{ii})
            freeparametercoordinatep1=length(sol.parameters)+1;
            sol.parameters=[sol.parameters jumpargument{ii}];
            OCMATINDIF.entrytimecoordinate{ii}=freeparametercoordinatep1:length(sol.parameters);
            if variationalcalculation
                %                 if length(jumpid)~=length(vjumpargument)
                %                     ocmaterror('Number of variational jump arguments and identifiers are unequal.')
                %                 end

                % initialize with zero
                freeparametercoordinatep1=length(sol.parameters)+1;
                sol.parameters=[sol.parameters vjumpargument{ii}];
                OCMATINDIF.ventrytimecoordinate{ii}=freeparametercoordinatep1:length(sol.parameters);
            end
        end
        if entryindex{ii}(end)==length(arcintv)
            OCMATINDIF.transversalityconditioncs{ii}=1;
        else
            OCMATINDIF.transversalityconditioncs{ii}=0;
        end
    elseif any([OCMATINDIF.stateconstraint{:}])
        OCMATINDIF.jumpargument{ii}=[];
        OCMATINDIF.entryindex{ii}=[];
        OCMATINDIF.jumpid{ii}=[];
        OCMATINDIF.transversalityconditioncs{ii}=0;
    else
        OCMATINDIF.transversalityconditioncs{ii}=0;
    end
end
depvar=dependentvar(ocTrjMP(1));
if ~isempty(freetimevector)
    arcintv=arcinterval(ocTrjMP(1));
    sol.parameters=[sol.parameters arcintv(end)];
    OCMATINDIF.freetimevectorcoordinate=length(sol.parameters);
end
if ~isempty(freeparametervector)
    if length(OCMATINDIF.freeparameterindex)==1 || isempty(OCMATINDIF.freeparametertarget)
        freeparametercoordinatep1=length(sol.parameters)+1;
        sol.parameters=[sol.parameters OCMATINDIF.parametervalue(OCMATINDIF.freeparameterindex)];
        OCMATINDIF.freeparametervectorcoordinate=freeparametercoordinatep1:length(sol.parameters);
        OCMATINDIF.freeinitialparametervalue=[];
    else
        OCMATINDIF.freeinitialparametervalue=OCMATINDIF.parametervalue(OCMATINDIF.freeparameterindex);
        OCMATINDIF.freeparameterdirection=OCMATINDIF.freeparametertarget-OCMATINDIF.freeinitialparametervalue;
        sol.parameters=[sol.parameters 0];
        OCMATINDIF.freeparametervectorcoordinate=length(sol.parameters);
    end
end
if static
    OCMATINDIF.staticparameterindex=staticparameterindex;
    for ii=1:indifforder
        freeparametercoordinatep1=length(sol.parameters)+1;
        par=modelparameter(ocTrjMP(ii));
        sol.parameters=[sol.parameters par(staticparameterindex)];
        OCMATINDIF.staticparametercoordinate{ii}=freeparametercoordinatep1:length(sol.parameters);
    end
end
if ~isempty(freestatevector)
    freestatevectorcoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters zeros(1,size(freestatevector,2))];
    OCMATINDIF.freestatevectorcoordinate=freestatevectorcoord:length(sol.parameters);
end
if ~isempty(fixinitstate)
    OCMATINDIF.initstate=depvar(fixinitstate,1);
else
    OCMATINDIF.initstate=[];    
end

if OCMATINDIF.continuationtype==0 % initial state continuation
    sol.parameters=[sol.parameters 0];
    OCMATINDIF.continuationcoordinate=length(sol.parameters);
    OCMATINDIF.continuationvector=OCMATINDIF.continuationtarget-OCMATINDIF.initialstate(continuationindex);
    OCMATINDIF.endtime=arcintv(end);
    OCMATINDIF.continuationindex=continuationindex;
elseif OCMATINDIF.continuationtype==1 % arctime continuation
    sol.parameters=[sol.parameters 0]; 
    OCMATINDIF.continuationcoordinate=length(sol.parameters);
    OCMATINDIF.continuationvector=OCMATINDIF.continuationtarget-arcintv(end);
    OCMATINDIF.continuationindex=continuationindex;
    OCMATINDIF.initialendtime=arcintv(end);
elseif OCMATINDIF.continuationtype==2
    OCMATINDIF.initialparameter=OCMATINDIF.parametervalue(parindex);
    sol.parameters=[sol.parameters 0];
    OCMATINDIF.continuationcoordinate=length(sol.parameters);
    OCMATINDIF.continuationvector=OCMATINDIF.continuationtarget.'-OCMATINDIF.initialparameter;
    OCMATINDIF.continuationindex=parindex;
    OCMATINDIF.endtime=arcintv(end);
end
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
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;
if variationalcalculation
    OCMATINDIF.variationaldynamicscoordinate=numberofodes+1:numberofodes+OCMATINDIF.variationalnumberofstates;
    OCMATINDIF.variationalstatecoordinate=dimensioncanonicalsystem+scoord(:).';
    OCMATINDIF.variationalcostatecoordinate=dimensioncanonicalsystem+cscoord(:).';

    numberofodes=numberofodes+dimensioncanonicalsystem;
    if size(sol.y,1)<numberofodes
        sol.y(OCMATINDIF.variationaldynamicscoordinate,:)=zeros(dimensioncanonicalsystem,length(sol.x));
    end
else
    if size(sol.y,1)>=2*dimensioncanonicalsystem
        sol.y(dimensioncanonicalsystem+1:2*dimensioncanonicalsystem,:)=[];
    end
end
OCMATINDIF.objectivevaluecoord=[];
if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATINDIF.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
end
OCMATINDIF.variationalobjectivevaluecoord=[];
if objectivevaluecalc & variationalcalculation
    numberofodes=numberofodes+1;
    OCMATINDIF.variationalobjectivevaluecoord=numberofodes;
    %sol.y(OCMATINDIF.variationalobjectivevaluecoord,:)=0;
end

if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATINDIF.exogenousnumberofstates;
else
    OCMATINDIF.exogenousdynamicscoordinate=[];
end

% F canonical system, FV variational system, FO objective function, FE
% exogenous functions, X states/costates, V derivative states/costates, E
% exogenous variables
OCMATINDIF.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATINDIF.dFDV=[]; % derivative of the canonical system with respect to the variational variable
OCMATINDIF.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATINDIF.dFVDX=[]; % derivative of the variatonal dynamics with respect to the dynamic variable
OCMATINDIF.dFVDO=[]; % derivative of the variatonal dynamics with respect to the objective variable
OCMATINDIF.dFVDE=[]; % derivative of the variatonal dynamics with respect to the exogenous variable
OCMATINDIF.dFVDV=[]; % derivative of the variatonal dynamics with respect to the variatonal variable
OCMATINDIF.dFVDPAR=[]; % derivative of the variatonal dynamics with respect to the parameters
OCMATINDIF.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODV=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATINDIF.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATINDIF.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDV=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATINDIF.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
OCMATINDIF.dFDX=zeros(dimensioncanonicalsystem);
OCMATINDIF.dFDXcoord1=1:dimensioncanonicalsystem;
OCMATINDIF.dFDXcoord2=1:dimensioncanonicalsystem;
coord1=dimensioncanonicalsystem;

% coord1 ... coordinates of canonical system / variational system /
% objective dynamics / exogenous dynamics
coord2=0;
if variationalcalculation
    % 1:n ... states, n+1:2n ... costates, 2n+1:3n .... derivative states,
    % 3n+1:4n ... derivative costates
    % 4n+1 ... objective value
    % 4n+2:4n+1+r ... exogenous values
    OCMATINDIF.dFDV=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDX=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDXcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATINDIF.dFVDXcoord2=coord2+1:2*dimensioncanonicalsystem;
    OCMATINDIF.dFVDV=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDVcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATINDIF.dFVDVcoord2=coord2+1:coord2+dimensioncanonicalsystem;
    coord1=coord1+dimensioncanonicalsystem;
    if exogenousfunction
        OCMATINDIF.dFVDE=zeros(dimensioncanonicalsystem,OCMATINDIF.exogenousnumberofstates);
    end
    if OCMATINDIF.objectivevaluecalc
        OCMATINDIF.dFVDO=zeros(dimensioncanonicalsystem,1);
    end
    OCMATINDIF.dFVDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
end
if objectivevaluecalc
    OCMATINDIF.dFDO=zeros(dimensioncanonicalsystem,1);
    OCMATINDIF.dFODO=0;
    OCMATINDIF.dFODX=zeros(1,dimensioncanonicalsystem);
    OCMATINDIF.dFODXcoord1=coord1+1:coord1+1;
    OCMATINDIF.dFODXcoord2=1:dimensioncanonicalsystem;
    coord1=coord1+1;
    if variationalcalculation
        OCMATINDIF.dFODV=zeros(1,dimensioncanonicalsystem);
    end
    if exogenousfunction
        OCMATINDIF.dFODE=zeros(1,OCMATINDIF.exogenousnumberofstates);
    end
    OCMATINDIF.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATINDIF.exogenousfunction
    OCMATINDIF.dFDE=zeros(dimensioncanonicalsystem,OCMATINDIF.exogenousnumberofstates);
    OCMATINDIF.dFEDE=zeros(OCMATINDIF.exogenousnumberofstates);
    OCMATINDIF.dFEDX=zeros(OCMATINDIF.exogenousnumberofstates,dimensioncanonicalsystem);
    OCMATINDIF.dFEDXcoord1=coord1+1:coord1+OCMATINDIF.exogenousnumberofstates;
    OCMATINDIF.dFEDXcoord2=1:dimensioncanonicalsystem;
    if variationalcalculation
        OCMATINDIF.dFEDV=zeros(OCMATINDIF.exogenousnumberofstates,dimensioncanonicalsystem);
    end
    if OCMATINDIF.objectivevaluecalc
        OCMATINDIF.dFEDO=zeros(OCMATINDIF.exogenousnumberofstates,1);
    end
    OCMATINDIF.dFEDPAR=zeros(OCMATINDIF.exogenousnumberofstates,length(sol.parameters));
end
OCMATINDIF.Jpar=[OCMATINDIF.dFDPAR;OCMATINDIF.dFVDPAR;OCMATINDIF.dFODPAR;OCMATINDIF.dFEDPAR];
OCMATINDIF.JX=[[OCMATINDIF.dFDX OCMATINDIF.dFDV OCMATINDIF.dFDO OCMATINDIF.dFDE];[OCMATINDIF.dFVDX OCMATINDIF.dFVDV OCMATINDIF.dFVDO OCMATINDIF.dFVDE];[OCMATINDIF.dFODX OCMATINDIF.dFODV OCMATINDIF.dFODO OCMATINDIF.dFODE];[OCMATINDIF.dFEDX OCMATINDIF.dFEDV OCMATINDIF.dFEDO OCMATINDIF.dFEDE]];

OCMATINDIF.ODEcoord=1:size(OCMATINDIF.Jpar,1);


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
