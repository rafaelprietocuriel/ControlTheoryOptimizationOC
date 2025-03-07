function sol=initocmat_FTE_P(ocObj,ocTrj,parindex,initialcoordinate,varargin)
%
% INITOCMAT_FTE_P initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_FTE_P(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see initocmat_FTE_T) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
fixendstate=[];
fixinitstate=[];
fixswitchingtimeindex=[];
followconstraint=[];
followobjectivevalue=[];
targetparametervalue=[];
optimalhorizon=false;
maxhorizon=[];
objectivevaluecalc=[];
hitvalue=[];
findoptimalparameter=[];
userbc=[];
exogenousfunction=[];
jumpargument=[];
jumpid=[];
entryindex=[];
freeparameter=[];
dynamicderivative=[]; % calculate the derivative of the solution path with respect to the specified continuation parameter
variationalcalculation=[];
vjumpargument=[];
vfreetime=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('oc trajectory is empty.')
    return
end
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixswitchingtimeindexidx=find(strcmpi(varargin,'fixswitchingtimeindex'));
followconstraintidx=find(strcmpi(varargin,'followconstraint'));
followobjectivevalueidx=find(strcmpi(varargin,'followobjectivevalue'));
userbcidx=find(strcmpi(varargin,'userbc'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
maxhorizonidx=find(strcmpi(varargin,'maxhorizon'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
findoptimalparameteridx=find(strcmpi(varargin,'findoptimalparameter'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
jumpididx=find(strcmpi(varargin,'jumpid'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
dynamicderivativeidx=find(strcmpi(varargin,'dynamicderivative'));
variationalcalculationidx=find(strcmpi(varargin,'variationalcalculation'));
vjumpargumentidx=find(strcmpi(varargin,'variationaljumpargument'));
vfreetimeidx=find(strcmpi(varargin,'variationalfreetime'));

if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end
if ~isempty(fixswitchingtimeindexidx)
    fixswitchingtimeindex=varargin{fixswitchingtimeindexidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(followconstraintidx)
    followconstraint=varargin{followconstraintidx+1};
end
if ~isempty(followobjectivevalueidx)
    followobjectivevalue=varargin{followobjectivevalueidx+1};
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
end
if ~isempty(findoptimalparameteridx)
    findoptimalparameter=varargin{findoptimalparameteridx+1};
end
if ~isempty(maxhorizonidx)
    maxhorizon=varargin{maxhorizonidx+1};
end
if isempty(maxhorizon)
    maxhorizon=inf;
end
if ~isempty(hitvalueidx)
    hitvalue=varargin{hitvalueidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
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
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(dynamicderivativeidx)
    dynamicderivative=varargin{dynamicderivativeidx+1};
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
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if ~isempty(followobjectivevalue)
    objectivevaluecalc=1;
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(parindex)
    ocmatmsg('Parameter is unknown.\n')
    return
end
if isempty(followconstraint)
    followconstraint=0;
end
if isempty(findoptimalparameter)
    findoptimalparameter=false;
end
if findoptimalparameter && (isempty(objectivevaluecalc) || objectivevaluecalc==0)
    objectivevaluecalc=true;
end
OCMATFTE.variationalcalculation=variationalcalculation;
if variationalcalculation
    OCMATFTE.variationalinitialstates=zeros(2*statenum(ocObj),1);
    OCMATFTE.variationalnumberofstates=length(OCMATFTE.variationalinitialstates);
end

OCMATFTE.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATFTE.exogenousinitialstates=exogenousinitialstates;
    OCMATFTE.exogenousnumberofstates=length(exogenousinitialstates);
end
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
if ~isempty(hitvalue)
    OCMATFTE.hitvaluefunc=funch{5}{4};
end
if ~isempty(userbc)
    OCMATFTE.userbcfunc=funch{5}{15};
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

if objectivevaluecalc
    OCMATFTE.objectivefunction=funch{8}{1};
    OCMATFTE.objectivefunctionjacobian=funch{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATFTE.objectivefunctionderivativetime=funch{8}{4};
    OCMATFTE.salvagevalue=funch{5}{6};
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
    OCMATFTE.variationalhamiltonian=funch{4}{6};
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
    end
end
if ~isautonomous(ocObj)
    OCMATFTE.canonicalsystemderivativetime=funch{2}{3};
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

OCMATFTE.findoptimalparameter=findoptimalparameter;
OCMATFTE.derivative=[];
OCMATFTE.userbc=userbc;
OCMATFTE.freeparameter=freeparameter;
if ~isempty(freeparameter)
    OCMATFTE.freeparameter=parameterindex(ocObj,freeparameter);
end

scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATFTE.statecostatecoord=[scoord(:).' cscoord(:).'];
sInfo=solverinfo(ocTrj);

sol=generatesolstruct(ocTrj,solver(ocTrj));

% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);
OCMATFTE.fixswitchingtimeindex=fixswitchingtimeindex;
if ~isempty(fixswitchingtimeindex)
    OCMATFTE.fixswitchingtime=sol.arcinterval(fixswitchingtimeindex);
end
OCMATFTE.freeswitchingtimeindex=2:length(sol.arcinterval(1:end-1));
OCMATFTE.freeswitchingtimeindex=setdiff(OCMATFTE.freeswitchingtimeindex,OCMATFTE.fixswitchingtimeindex);
numswitchtimes=length(sol.arcinterval(2:end))-1;
sol.parameters=sol.arcinterval(OCMATFTE.freeswitchingtimeindex);
freeparametercoordinatep1=length(sol.parameters)+1;
OCMATFTE.initarcinterval=sol.arcinterval;
OCMATFTE.switchingtimetype=zeros(1,length(sol.arcinterval));
OCMATFTE.switchingtimetype(OCMATFTE.freeswitchingtimeindex)=1;
OCMATFTE.freeswitchingtimecoordinate=1:length(sol.parameters);
if variationalcalculation
    if length(vfreetime)~=length(OCMATFTE.freeswitchingtimeindex)
        ocmaterror('Number of variational free time arguments and number of free time arguments are different.')
    end
    sol.parameters=[sol.parameters vfreetime];
    OCMATFTE.vfreetimecoord=freeparametercoordinatep1:length(sol.parameters);
end
if optimalhorizon
    sol.parameters=[sol.parameters sol.arcinterval(end)];
    OCMATFTE.optimalhorizoncoord=length(sol.parameters);
    if variationalcalculation
        if isfield(sInfo,'voptimalhorizoncoord')
            % use from previous solution
            sol.parameters=[sol.parameters sInfo.parameter(sInfo.voptimalhorizoncoord)];
        else
            % initialize with zero
            sol.parameters=[sol.parameters 0];
        end
        OCMATFTE.voptimalhorizoncoord=length(sol.parameters);
    end
else
    OCMATFTE.truncationtime=sol.arcinterval(end);
end
OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.maxhorizon=maxhorizon;
OCMATFTE.constraintidx=followconstraint;
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
            freeparametercoordinatep1=length(sol.parameters)+1;
        end
    end
    if entryindex(end)==length(sol.arcinterval)
        OCMATFTE.transversalityconditioncs=1;
    end
end
if ~isempty(freeparameter)
    freeparametercoordinate=length(sol.parameters)+1;
    sol.parameters=[sol.parameters OCMATFTE.parametervalue(OCMATFTE.freeparameter)];
    OCMATFTE.freeparametercoordinate=freeparametercoordinate:length(sol.parameters);
%     if variationalcalculation
%         if isfield(sInfo,'vfreeparametercoordinate')
%             % use from previous solution
%             sol.parameters=[sol.parameters sInfo.parameter(sInfo.vfreeparametercoordinate)];
%         else
%             % initialize with zero
%             sol.parameters=[sol.parameters zeros(1,length(OCMATFTE.freeparameter))];
%         end
%         OCMATFTE.vfreeparametercoordinate=freeparametercoordinatep1:length(sol.parameters);
%     end
end

OCMATFTE.continuationcoordinate=length(sol.parameters)+1;
OCMATFTE.initialparametervalue=OCMATFTE.parametervalue(parindex);
if isempty(targetparametervalue)
    OCMATFTE.continuationvector=OCMATFTE.parametervalue(parindex);
else
    if length(parindex)==length(targetparametervalue)
        OCMATFTE.continuationvector=targetparametervalue(:).'-OCMATFTE.parametervalue(parindex);
    else
        ocmatmsg('Free parameter values and target values do not agree.\n')
        return
    end
end
sol.parameters=[sol.parameters 0];

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
% if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj) || objectivevaluecalc>1
%     OT=discountedsalvagevalue(ocObj,ocTrj);
%     o=objectivefunction(ocObj,ocTrj,1);
%     sol.y(end+1,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
% end
% if dynamicderivative
%     J=jacobian(ocObj,ocTrj,[],[],1);
%     if objectivevaluecalc 
%         JO=objectivefunctionjacobian(ocObj,ocTrj,1);
%     end
%     if objectivevaluecalc 
%         JE=exogenousfunctionjacobian(ocObj,ocTrj,1);
%     end
%     
% end
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
if isempty(followconstraint) && isempty(fixinitstate)
    fixinitstate=initialcoordinate;
end

OCMATFTE.followobjectivevalue=followobjectivevalue; % value of the objective to be followd

OCMATFTE.initialtime=sol.x0;
OCMATFTE.initialstate=depvar(initialcoordinate,1);
OCMATFTE.switchtimecoord=1:numswitchtimes;

OCMATFTE.continuationindex=parindex;
OCMATFTE.targetparametervalue=targetparametervalue;
%OCMATFTE.targetvalue=targetparametervalue;

OCMATFTE.objectivevaluecalc=objectivevaluecalc;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitstatecoord=fixinitstate;
if ~isempty(fixendstate)
    OCMATFTE.endstate=depvar(fixendstate,end);
end
if ~isempty(fixinitstate)
    OCMATFTE.initstate=depvar(fixinitstate,1);
else
    OCMATFTE.initstate=[];    
end
pathname=OCMATFTE.datapath();
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;
if variationalcalculation
    OCMATFTE.variationaldynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.variationalnumberofstates;
    OCMATFTE.variationalstatecoordinate=dimensioncanonicalsystem+scoord(:).';
    OCMATFTE.variationalcostatecoordinate=dimensioncanonicalsystem+cscoord(:).';

    numberofodes=numberofodes+dimensioncanonicalsystem;
    if size(sol.y,1)<numberofodes
        sol.y(OCMATFTE.variationaldynamicscoordinate,:)=zeros(dimensioncanonicalsystem,length(sol.x));
    end
else
    if size(sol.y,1)>=2*dimensioncanonicalsystem
        sol.y(dimensioncanonicalsystem+1:2*dimensioncanonicalsystem,:)=[];
    end
end
OCMATFTE.objectivevaluecoord=[];
if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATFTE.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
end
if objectivevaluecalc && length(sol.y(:,1))<numberofodes || objectivevaluecalc>1
    OT=discountedsalvagevalue(ocObj,ocTrj);
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(OCMATFTE.objectivevaluecoord,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.exogenousnumberofstates;
    numberofodes=numberofodes+OCMATFTE.exogenousnumberofstates;
else
    OCMATFTE.exogenousdynamicscoordinate=[];
end

OCMATFTE.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATFTE.dFDV=[]; % derivative of the canonical system with respect to the variational variable
OCMATFTE.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
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
OCMATFTE.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATFTE.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
OCMATFTE.dFDX=zeros(dimensioncanonicalsystem);
OCMATFTE.dFDXcoord1=1:dimensioncanonicalsystem;
OCMATFTE.dFDXcoord2=1:dimensioncanonicalsystem;
coord1=dimensioncanonicalsystem;
coord2=0;
if variationalcalculation
    OCMATFTE.dFDV=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDX=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDXcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATFTE.dFVDXcoord2=coord2+1:2*dimensioncanonicalsystem;
    OCMATFTE.dFVDV=zeros(dimensioncanonicalsystem);
    OCMATFTE.dFVDVcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATFTE.dFVDVcoord2=coord2+1:coord2+dimensioncanonicalsystem;
    coord1=coord1+dimensioncanonicalsystem;
    if exogenousfunction
        OCMATFTE.dFVDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    end
    if OCMATFTE.objectivevaluecalc
        OCMATFTE.dFVDO=zeros(dimensioncanonicalsystem,1);
    end
    OCMATFTE.dFVDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
end
if objectivevaluecalc
    OCMATFTE.dFDO=zeros(dimensioncanonicalsystem,1);
    OCMATFTE.dFODO=0;
    OCMATFTE.dFODX=zeros(1,dimensioncanonicalsystem);
    OCMATFTE.dFODXcoord1=coord1+1:coord1+1;
    OCMATFTE.dFODXcoord2=1:dimensioncanonicalsystem;
    coord1=coord1+1;
    if variationalcalculation
        OCMATFTE.dFODV=zeros(1,dimensioncanonicalsystem);
    end
    if exogenousfunction
        OCMATFTE.dFODE=zeros(1,OCMATFTE.exogenousnumberofstates);
    end
    OCMATFTE.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATFTE.exogenousfunction
    OCMATFTE.dFDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDE=zeros(OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDX=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    OCMATFTE.dFEDXcoord1=coord1+1:coord1+OCMATFTE.exogenousnumberofstates;
    OCMATFTE.dFEDXcoord2=1:dimensioncanonicalsystem;
    if variationalcalculation
        OCMATFTE.dFEDV=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    end
    if OCMATFTE.objectivevaluecalc
        OCMATFTE.dFEDO=zeros(OCMATFTE.exogenousnumberofstates,1);
    end
    OCMATFTE.dFEDPAR=zeros(OCMATFTE.exogenousnumberofstates,length(sol.parameters));
end
OCMATFTE.Jpar=[OCMATFTE.dFDPAR;OCMATFTE.dFVDPAR;OCMATFTE.dFODPAR;OCMATFTE.dFEDPAR];
 
OCMATFTE.Jext=[[OCMATFTE.dFDV OCMATFTE.dFDO OCMATFTE.dFDE];[OCMATFTE.dFVDV OCMATFTE.dFVDO OCMATFTE.dFVDE];[OCMATFTE.dFODV OCMATFTE.dFODO OCMATFTE.dFODE];[OCMATFTE.dFEDV OCMATFTE.dFEDO OCMATFTE.dFEDE]];
OCMATFTE.JX=[OCMATFTE.dFDX;OCMATFTE.dFVDX;OCMATFTE.dFODX;OCMATFTE.dFEDX];
OCMATFTE.ODEcoord=1:size(OCMATFTE.Jpar,1);

[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);

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