function sol=initocmat_FTE(ocObj,ocTrj,contcoordinate,targetvalue,varargin)
%
% if ismultimodel(ocObj)
%     sol=multimodel/initocmat_FTE(ocObj,ocTrj,contcoordinate,targetvalue,varargin{:});
%     return
% end
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
fixendstate=[];
freeendstate=[];
fixinitstate=[];
freeinitstate=[];
initstateinequality=[]; % value of the Lagrangian
freeparameter=0;
initialhorizon=[];
maxhorizon=[];
objectivevaluecalc=[];
hitvalue=[];
userbc=[];
jumpargument=[];
jumpid=[];
entryindex=[];
exogenousfunction=[];
variationalcalculation=[];
includevariationalobjectivevalue=[];
freevector=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('oc trajectory is empty.')
    return
end
maxhorizonidx=find(strcmpi(varargin,'maxhorizon'));
optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
initialhorizonidx=find(strcmpi(varargin,'initialhorizon'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
freeendstateidx=find(strcmpi(varargin,'freeendstate'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
freeinitstateidx=find(strcmpi(varargin,'freeinitstate'));
initstateinequalityidx=find(strcmpi(varargin,'initstateinequality'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
userbcidx=find(strcmpi(varargin,'userbc'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
jumpididx=find(strcmpi(varargin,'jumpid'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
variationalcalculationidx=find(strcmpi(varargin,'variationalcalculation'));
vjumpargumentidx=find(strcmpi(varargin,'variationaljumpargument'));
vfreetimeidx=find(strcmpi(varargin,'variationalfreetime'));
includevariationalobjectivevalueidx=find(strcmpi(varargin,'includevariationalobjectivevalue'));
freevectoridx=find(strcmpi(varargin,'freevector'));

if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
else
    optimalhorizon=0;
end
if ~isempty(maxhorizonidx)
    maxhorizon=varargin{maxhorizonidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(freeendstateidx)
    freeendstate=varargin{freeendstateidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end
if ~isempty(freeinitstateidx)
    freeinitstate=varargin{freeinitstateidx+1};
end
if ~isempty(initstateinequalityidx)
    initstateinequality=varargin{initstateinequalityidx+1};
end
if ~isempty(initialhorizon)
    initialhorizon=varargin{initialhorizonidx+1};
else
    arcinterval(ocTrj);
    initialhorizon=ans(end);
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
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
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if ~isempty(variationalcalculationidx)
    variationalcalculation=varargin{variationalcalculationidx+1};
end
if ~isempty(includevariationalobjectivevalueidx)
    includevariationalobjectivevalue=varargin{includevariationalobjectivevalueidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end

if ~isempty(vjumpargumentidx)
    vjumpargument=varargin{vjumpargumentidx+1};
end
if ~isempty(vfreetimeidx)
    vfreetime=varargin{vfreetimeidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(variationalcalculation)
    variationalcalculation=false;
end
if isempty(includevariationalobjectivevalue)
    includevariationalobjectivevalue=false;
end

if isempty(maxhorizon)
    maxhorizon=inf;
end
targetvalue=targetvalue(:);

if freeparameter
    parindex=varargin{freeparameteridx+1};
    if ischar(parindex)
        parindex=parameterindex(ocObj,parindex);
    end
    OCMATFTE.freeparameter=1;
    OCMATFTE.parameterindex=parindex;
else
    OCMATFTE.freeparameter=0;
end
OCMATFTE.userbc=userbc;
OCMATFTE.userbcvalue=userbc;

OCMATFTE.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATFTE.exogenousinitialstates=exogenousinitialstates;
    OCMATFTE.exogenousnumberofstates=length(exogenousinitialstates);
end
OCMATFTE.variationalcalculation=variationalcalculation;
if OCMATFTE.variationalcalculation
    OCMATFTE.variationalinitialstates=zeros(2*statenum(ocObj),1);
    OCMATFTE.variationalnumberofstates=length(OCMATFTE.variationalinitialstates);
end
% test if pure state constraints are defined and costates can jump (entry
% times)
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
if ~isempty(userbc)
    %OCMATFTE.userbc=funch{5}{8};
    OCMATFTE.userbcfunc=funch{5}{15};
end
if ~isempty(freeinitstate)
    OCMATFTE.bcfreeinitial=funch{5}{9};
end
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

OCMATFTE.includevariationalobjectivevalue=includevariationalobjectivevalue;
if objectivevaluecalc
    OCMATFTE.objectivefunction=funch{8}{1};
    OCMATFTE.objectivefunctionjacobian=funch{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATFTE.objectivefunctionderivativetime=funch{8}{4};
    OCMATFTE.salvagevalue=funch{5}{6};
    if OCMATFTE.includevariationalobjectivevalue
        OCMATFTE.variationalobjectivefunction=funch{8}{5};
        OCMATFTE.variationalobjectivefunctionjacobian=funch{8}{6};
        OCMATFTE.variationalobjectivefunctionparameterjacobian=funch{8}{7};
        OCMATFTE.variationalobjectivefunctionderivativetime=funch{8}{8};
        OCMATFTE.variationalsalvagevalue=funch{8}{9};
    end
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamics=funch{4}{1};
    OCMATFTE.exogenousjacobian=funch{4}{2};
    OCMATFTE.exogenousparameterjacobian=funch{4}{3};
end
if OCMATFTE.variationalcalculation
    OCMATFTE.variationaldynamics=funch{4}{4};
    OCMATFTE.variationaljacobian=funch{4}{5};
    OCMATFTE.variationalparameterjacobian=funch{4}{6};
    OCMATFTE.variationalhamiltonian=funch{4}{10};
    OCMATFTE.variationalguard=funch{5}{9};
    OCMATFTE.variationalreset=funch{5}{10};
    OCMATFTE.variationalbcinitial=funch{5}{11};
    OCMATFTE.variationalbctransversality=funch{5}{12};
    if OCMATFTE.stateconstraint
        OCMATFTE.variationalbcstateconstraint=funch{5}{13};
        OCMATFTE.variationalbctransversalitysc=funch{5}{14};
    end
    if OCMATFTE.exogenousfunction
        OCMATFTE.exogenousvariationaldynamics=funch{4}{7};
        OCMATFTE.exogenousjacobian4variationalargument=funch{4}{8};
        OCMATFTE.exogenousvariationaldynamicsjacobian=funch{4}{9};
    end
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

%sol=generatesolstruct(ocTrj,solver(ocTrj));
solold=generatesolstruct(ocTrj,solver(ocTrj));
sol=solold;
sol.y=[];
sol.x=solold.x;
sol.arcarg=solold.arcarg;
sol.arcinterval=solold.arcinterval;
sol.arcposition=solold.arcposition;
sol.x0=solold.x0;

OCMATFTE.freevector=freevector;

OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.maxhorizon=maxhorizon;
% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATFTE.statecostatecoord=[scoord(:).' cscoord(:).'];
OCMATFTE.costatecoord=cscoord(:).';

numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
end
OCMATFTE.freeswitchingtimeindex=2:length(sol.arcinterval(1:end-1));
if optimalhorizon
    sol.parameters=[sol.parameters initialhorizon];
end
OCMATFTE.switchtimecoordinate=1:numswitchtimes;
if  OCMATFTE.freeparameter
    %OCMATFTE.freevectorcoord=[];
    sol.parameters=[sol.parameters OCMATFTE.parametervalue(parindex)];
    OCMATFTE.parametercoord=length(sol.parameters)-length(parindex)+(1:length(parindex));
end
freeparametercoordinatep1=length(sol.parameters)+1;

if variationalcalculation
    if length(vfreetime)~=length(OCMATFTE.freeswitchingtimeindex)
        ocmaterror('Number of variational free time arguments and number of free time arguments are different.')
    end
    sol.parameters=[sol.parameters vfreetime];
    OCMATFTE.vfreetimecoord=freeparametercoordinatep1:length(sol.parameters);
end

OCMATFTE.transversalityconditioncs=0;
if OCMATFTE.stateconstraint
    OCMATFTE.jumpid=zeros(1,length(sol.arcinterval));
    OCMATFTE.jumpargument=jumpargument;
    if length(jumpid)~=length(jumpargument)
        ocmaterror('Number of jump arguments and identifiers are unequal.')
    end
    OCMATFTE.entryindex=zeros(1,length(sol.arcinterval));
    OCMATFTE.entryindex(entryindex)=1:length(jumpargument);
    OCMATFTE.jumpid(entryindex)=jumpid;
    if ~isempty(jumpargument)
        entrytimecoordinate=length(sol.parameters)+1;
        sol.parameters=[sol.parameters jumpargument];
        OCMATFTE.entrytimecoordinate=entrytimecoordinate:length(sol.parameters);
        freeparametercoordinatep1=length(sol.parameters)+1;
        if variationalcalculation
            if length(jumpid)~=length(vjumpargument)
                ocmaterror('Number of variational jump arguments and identifiers are unequal.')
            end

            % initialize with zero
            sol.parameters=[sol.parameters vjumpargument];
            OCMATFTE.ventrytimecoordinate=freeparametercoordinatep1:length(sol.parameters);
        end
    end
    if entryindex(end)==length(sol.arcinterval)
        OCMATFTE.transversalityconditioncs=1;
    end
end
OCMATFTE.initstateinequality=initstateinequality;
OCMATFTE.initstateinequalitycoord=[];
if ~isempty(initstateinequality)
    initstateinequalitycoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters initstateinequality];
    OCMATFTE.initstateinequalitycoord=initstateinequalitycoord:length(sol.parameters);
end
if ~isempty(freevector)
    freevectorcoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters zeros(1,size(freevector,2))];
    OCMATFTE.freevectorcoord=freevectorcoord:length(sol.parameters);
else
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


OCMATFTE.initialtime=sol.x0;
if optimalhorizon
    if numswitchtimes
        if ~isempty(OCMATFTE.switchtimecoordinate)
        OCMATFTE.optimalhorizoncoord=OCMATFTE.switchtimecoordinate(end)+1;
        else
            OCMATFTE.optimalhorizoncoord=1;
        end
    else
        OCMATFTE.optimalhorizoncoord=1;
    end
else
    OCMATFTE.endtime=sol.arcinterval(end);
    OCMATFTE.optimalhorizoncoord=[];
end
OCMATFTE.targetcoordinate=contcoordinate;
%if isempty(fixinitstate)
%    fixinitstate=contcoordinate;
%end
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);

if isempty(fixinitstate) %statenum(ocObj)>=length(fixinitstate)
    if isempty(freeinitstate)
        if ~isempty(fixinitstate)
            freeinitstate=setdiff(statenum(ocObj),fixinitstate);
        else
            freeinitstate=[];
        end
    end
    if ~isempty(contcoordinate)
        OCMATFTE.startvalue=depvar(contcoordinate,1);
    else
        OCMATFTE.startvalue=depvar(1:dimensioncanonicalsystem/2,1);
        OCMATFTE.targetcoordinate=1:dimensioncanonicalsystem/2;
    end
    OCMATFTE.continitstate=true;
elseif isempty(fixendstateidx)
    OCMATFTE.startvalue=depvar(contcoordinate,end);
    OCMATFTE.continitstate=false;
end
if ~isempty(contcoordinate)
    OCMATFTE.continuationvector=targetvalue-OCMATFTE.startvalue;
end
OCMATFTE.objectivevaluecalc=objectivevaluecalc;
OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.fixinitstatecoord=fixinitstate;
OCMATFTE.fixendstatecoord=fixendstate;
if ~isempty(fixendstate)
    OCMATFTE.endstate=depvar(fixendstate,end);
    %OCMATFTE.initstate=[];
end
if ~isempty(freeendstate)
    OCMATFTE.freeendstatecoord=freeendstate;
else
    if ~isempty(fixendstate)
        OCMATFTE.freeendstatecoord=fixendstate;
    else
        OCMATFTE.freeendstatecoord=1:statenum(ocObj);
    end
end
OCMATFTE.freeinitstatecoord=freeinitstate;
if ~isempty(fixinitstate)
    OCMATFTE.initstate=depvar(fixinitstate,1);
    %OCMATFTE.endstate=[];
end
pathname=OCMATFTE.datapath();


% build up the generalized state variables
numberofodes=dimensioncanonicalsystem;
sol.y(1:dimensioncanonicalsystem,:)=solold.y(1:dimensioncanonicalsystem,:);

if ~isempty(hitvalue)
    %OCMATFTE.targetdistance=norm(OCMATFTE.startvalue-hitvalue)/norm(OCMATFTE.startvalue-targetvalue);
    OCMATFTE.hitvalue=hitvalue;
    try
        OCMATFTE.hitvaluefunc=funch{5}{4};
    catch
        OCMATFTE.hitvaluefunc=[];
    end
else
     OCMATFTE.targetdistance=[];
     OCMATFTE.hitvalue=[];
end

if variationalcalculation
    OCMATFTE.variationaldynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.variationalnumberofstates;
    OCMATFTE.variationalstatecoordinate=dimensioncanonicalsystem+scoord(:).';
    OCMATFTE.variationalcostatecoordinate=dimensioncanonicalsystem+cscoord(:).';

    numberofodes=numberofodes+dimensioncanonicalsystem;
    if size(solold.y,1)<numberofodes
        sol.y(OCMATFTE.variationaldynamicscoordinate,:)=zeros(dimensioncanonicalsystem,length(sol.x));
    else
        sol.y(OCMATFTE.variationaldynamicscoordinate,:)=solold.y(OCMATFTE.variationaldynamicscoordinate,:);
    end
end

OCMATFTE.objectivevaluecoord=[];

if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATFTE.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
end
if objectivevaluecalc && length(solold.y(:,1))<numberofodes
    OT=discountedsalvagevalue(ocObj,ocTrj);
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(OCMATFTE.objectivevaluecoord,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
else
    sol.y(OCMATFTE.objectivevaluecoord,:)=solold.y(OCMATFTE.objectivevaluecoord,:);
end
OCMATFTE.variationalobjectivevaluecoord=[];
if includevariationalobjectivevalue
    numberofodes=numberofodes+1;
    OCMATFTE.variationalobjectivevaluecoord=numberofodes;
    if size(solold.y,1)<numberofodes
        sol.y(OCMATFTE.variationalobjectivevaluecoord,:)=0;
    else
        sol.y(OCMATFTE.variationalobjectivevaluecoord,:)=solold.y(OCMATFTE.variationalobjectivevaluecoord,:);
    end
end

if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.exogenousnumberofstates;
else
    OCMATFTE.exogenousdynamicscoordinate=[];
end
% F canonical system, FV variational system, FO objective function, FE
% exogenous functions, X states/costates, V derivative states/costates, E
% exogenous variables
OCMATFTE.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATFTE.dFDV=[]; % derivative of the canonical system with respect to the variational variable
OCMATFTE.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATFTE.dFDPAR=[]; % derivative of the canonical system with respect to the parameters
OCMATFTE.dFVDX=[]; % derivative of the variatonal dynamics with respect to the dynamic variable
OCMATFTE.dFVDO=[]; % derivative of the variatonal dynamics with respect to the objective variable
OCMATFTE.dFVDE=[]; % derivative of the variatonal dynamics with respect to the exogenous variable
OCMATFTE.dFVDV=[]; % derivative of the variatonal dynamics with respect to the variatonal variable
OCMATFTE.dFVDPAR=[]; % derivative of the variatonal dynamics with respect to the parameters
OCMATFTE.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODV=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATFTE.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATFTE.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATFTE.dFEDV=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATFTE.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the parameters

OCMATFTE.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
OCMATFTE.dFDX=zeros(dimensioncanonicalsystem);
OCMATFTE.dFDXcoord1=1:dimensioncanonicalsystem;
OCMATFTE.dFDXcoord2=1:dimensioncanonicalsystem;
coord1=dimensioncanonicalsystem;

% coord1 ... coordinates of canonical system / variational system /
% objective dynamics / exogenous dynamics
coord2=0;
if variationalcalculation
    % 1:n ... states, n+1:2n ... costates, 2n+1:3n .... derivative states,
    % 3n+1:4n ... derivative costates
    % 4n+1 ... objective value
    % 4n+2:4n+1+r ... exogenous values
    OCMATFTE.dFDV=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDX=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDXcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    %OCMATFTE.dFVDXcoord2=coord2+1:2*dimensioncanonicalsystem;
    OCMATFTE.dFVDV=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDVcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATFTE.dFVDVcoord2=coord2+1:coord2+dimensioncanonicalsystem;
    coord1=coord1+dimensioncanonicalsystem;
    if exogenousfunction
        OCMATFTE.dFVDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    end
    if OCMATFTE.includevariationalobjectivevalue
        OCMATFTE.dFVDO=zeros(dimensioncanonicalsystem,2);
    elseif objectivevaluecalc
        OCMATFTE.dFVDO=zeros(dimensioncanonicalsystem,1);
    end
    OCMATFTE.dFVDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
end
if objectivevaluecalc
    if includevariationalobjectivevalue
        numo=2;
    else
        numo=1;
    end
    OCMATFTE.dFDO=zeros(dimensioncanonicalsystem,numo);
    OCMATFTE.dFODO=zeros(numo);
    OCMATFTE.dFODX=zeros(numo,dimensioncanonicalsystem);
    OCMATFTE.dFODXcoord1=coord1+1:coord1+numo;
    coord1=coord1+numo;
    OCMATFTE.dFODXcoord2=1:dimensioncanonicalsystem;
    if exogenousfunction
        OCMATFTE.dFODE=zeros(numo,OCMATFTE.exogenousnumberofstates);
    end
    if variationalcalculation
        OCMATFTE.dFODV=zeros(numo,dimensioncanonicalsystem);
    end
    OCMATFTE.dFODPAR=zeros(numo,length(sol.parameters));
end
if OCMATFTE.exogenousfunction
    OCMATFTE.dFDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDE=zeros(OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDX=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    OCMATFTE.dFEDXcoord1=coord1+1:coord1+OCMATFTE.exogenousnumberofstates;
    OCMATFTE.dFEDXcoord2=1:dimensioncanonicalsystem;
    coord1=coord1+OCMATFTE.exogenousnumberofstates;
    if variationalcalculation
        OCMATFTE.dFEDV=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    end
    if OCMATFTE.objectivevaluecalc
        OCMATFTE.dFEDO=zeros(OCMATFTE.exogenousnumberofstates,numo);
    end
    OCMATFTE.dFEDPAR=zeros(OCMATFTE.exogenousnumberofstates,length(sol.parameters));
end
OCMATFTE.Jpar=[OCMATFTE.dFDPAR;OCMATFTE.dFVDPAR;OCMATFTE.dFODPAR;OCMATFTE.dFEDPAR];
OCMATFTE.JX=[[OCMATFTE.dFDX OCMATFTE.dFDV OCMATFTE.dFDO OCMATFTE.dFDE];[OCMATFTE.dFVDX OCMATFTE.dFVDV OCMATFTE.dFVDO OCMATFTE.dFVDE];[OCMATFTE.dFODX OCMATFTE.dFODV OCMATFTE.dFODO OCMATFTE.dFODE];[OCMATFTE.dFEDX OCMATFTE.dFEDV OCMATFTE.dFEDO OCMATFTE.dFEDE]];

% OCMATFTE.dFDO=[]; % derivative of the canonical system with respect to the objective variable
% OCMATFTE.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
% OCMATFTE.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
% OCMATFTE.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
% OCMATFTE.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
% OCMATFTE.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
% OCMATFTE.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
% OCMATFTE.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
% OCMATFTE.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
% OCMATFTE.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
% 
% OCMATFTE.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
% OCMATFTE.dFDX=zeros(dimensioncanonicalsystem);
% OCMATFTE.dFDXcoord1=1:dimensioncanonicalsystem;
% OCMATFTE.dFDXcoord2=1:dimensioncanonicalsystem;
% if objectivevaluecalc
%     OCMATFTE.dFDO=zeros(dimensioncanonicalsystem,1);
%     OCMATFTE.dFODO=0;
%     OCMATFTE.dFODX=zeros(1,dimensioncanonicalsystem);
%     if exogenousfunction
%         OCMATFTE.dFODE=zeros(1,OCMATFTE.exogenousnumberofstates);
%     end
%     OCMATFTE.dFODPAR=zeros(1,length(sol.parameters));
% end
% if OCMATFTE.exogenousfunction
%     OCMATFTE.dFDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
%     OCMATFTE.dFEDE=zeros(OCMATFTE.exogenousnumberofstates);
%     OCMATFTE.dFEDX=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
%     if OCMATFTE.objectivevaluecalc
%         OCMATFTE.dFEDO=zeros(OCMATFTE.exogenousnumberofstates,1);
%     end
%     OCMATFTE.dFEDPAR=zeros(OCMATFTE.exogenousnumberofstates,length(sol.parameters));
% end
% OCMATFTE.Jpar=[OCMATFTE.dFDPAR;OCMATFTE.dFODPAR;OCMATFTE.dFEDPAR];
OCMATFTE.ODEcoord=1:size(OCMATFTE.Jpar,1);
 
%OCMATFTE.Jext=[[OCMATFTE.dFDO OCMATFTE.dFDE];[OCMATFTE.dFODO OCMATFTE.dFODE];[OCMATFTE.dFEDO OCMATFTE.dFEDE]];

[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);

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