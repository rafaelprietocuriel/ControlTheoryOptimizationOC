function sol=initocmat_FTE_J(ocObj,ocTrj,stateconstraintindex,contexitindex,contentryindex,varargin)
%
% if ismultimodel(ocObj)
%     sol=multimodel/initocmat_FTE(ocObj,ocTrj,contcoordinate,targetvalue,varargin{:});
%     return
% end
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
fixendstate=[];
fixinitstate=[];
freeparameter=0;
initialhorizon=[];
maxhorizon=[];
objectivevaluecalc=[];
hitvalue=[];
userbc=[];
jumpargument=[];
entryindex=[];
exogenousfunction=[];
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
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
hitvalueidx=find(strcmpi(varargin,'hitvalue'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
userbcidx=find(strcmpi(varargin,'userbc'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
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
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
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
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(maxhorizon)
    maxhorizon=inf;
end


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
    OCMATFTE.userbc=funch{5}{8};
end
if OCMATFTE.stateconstraint
    OCMATFTE.bcstateconstraint=funch{5}{7};
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
OCMATFTE.guard4jump=funch{7}{6};
OCMATFTE.jacobianguard=funch{7}{7};
OCMATFTE.jacobianreset=funch{7}{8};
OCMATFTE.domaindiscretization=funch{7}{9};

if objectivevaluecalc
    OCMATFTE.objectivefunction=funch{8}{1};
    OCMATFTE.objectivefunctionjacobian=funch{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATFTE.salvagevalue=funch{5}{6};
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamics=funch{4}{1};
    OCMATFTE.exogenousjacobian=funch{4}{2};
    OCMATFTE.exogenousparameterjacobian=funch{4}{3};
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

sol=generatesolstruct(ocTrj,solver(ocTrj));

OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.maxhorizon=maxhorizon;
% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATFTE.statecostatecoord=[scoord(:).' cscoord(:).'];

numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
end
if optimalhorizon
    sol.parameters=[sol.parameters initialhorizon];
end
OCMATFTE.switchtimecoord=1:numswitchtimes;
if  OCMATFTE.freeparameter
    %OCMATFTE.freevectorcoord=[];
    sol.parameters=[sol.parameters OCMATFTE.parametervalue(parindex)];
    OCMATFTE.parametercoord=length(sol.parameters)-length(parindex)+(1:length(parindex));
end
if OCMATFTE.stateconstraint
    OCMATFTE.jumpargument=jumpargument;
    OCMATFTE.entryindex=zeros(1,length(sol.arcinterval));
    entryindex=[entryindex contentryindex];
    OCMATFTE.entryindex(entryindex)=1:length(jumpargument);
    if ~isempty(jumpargument)
        entrytimecoordinate=length(sol.parameters)+1;
        sol.parameters=[sol.parameters jumpargument];
        OCMATFTE.entrytimecoordinate=entrytimecoordinate:length(sol.parameters);
    end
    l=lagrangemultiplier(ocObj,ocTrj,1);
    arcpos=arcposition(ocTrj);
    OCMATFTE.lagrangemultiplieratexit=l(stateconstraintindex,arcpos(2,contexitindex));
    OCMATFTE.exitindex=contexitindex;
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
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj) || objectivevaluecalc>1
    OT=discountedsalvagevalue(ocObj,ocTrj);
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(end+1,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
end


OCMATFTE.initialtime=sol.x0;
if optimalhorizon
    if numswitchtimes
        if ~isempty(OCMATFTE.switchtimecoord)
        OCMATFTE.optimalhorizoncoord=OCMATFTE.switchtimecoord(end)+1;
        else
            OCMATFTE.optimalhorizoncoord=1;
        end
    else
        OCMATFTE.optimalhorizoncoord=1;
    end
else
    OCMATFTE.truncationtime=sol.arcinterval(end);
    OCMATFTE.optimalhorizoncoord=[];
end
OCMATFTE.initstate=depvar(fixinitstate,1);
OCMATFTE.endstate=depvar(fixendstate,end);
OCMATFTE.objectivevaluecalc=objectivevaluecalc;
OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitstatecoord=fixinitstate;
pathname=OCMATFTE.datapath();
if objectivevaluecalc
    OCMATFTE.objectivevaluecoord=size(sol.y,1);
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(ii).numode,:);
end
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

dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;
OCMATFTE.objectivevaluecoord=[];
if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATFTE.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATFTE.exogenousnumberofstates;
else
    OCMATFTE.exogenousdynamicscoordinate=[];
end

OCMATFTE.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATFTE.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATFTE.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATFTE.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATFTE.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATFTE.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATFTE.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATFTE.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
if objectivevaluecalc
    OCMATFTE.dFDO=zeros(dimensioncanonicalsystem,1);
    OCMATFTE.dFODO=0;
    OCMATFTE.dFODX=zeros(1,dimensioncanonicalsystem);
    if exogenousfunction
        OCMATFTE.dFODE=zeros(1,OCMATFTE.exogenousnumberofstates);
    end
    OCMATFTE.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATFTE.exogenousfunction
    OCMATFTE.dFDE=zeros(dimensioncanonicalsystem,OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDE=zeros(OCMATFTE.exogenousnumberofstates);
    OCMATFTE.dFEDX=zeros(OCMATFTE.exogenousnumberofstates,dimensioncanonicalsystem);
    if OCMATFTE.objectivevaluecalc
        OCMATFTE.dFEDO=zeros(OCMATFTE.exogenousnumberofstates,1);
    end
    OCMATFTE.dFEDPAR=zeros(OCMATFTE.exogenousnumberofstates,length(sol.parameters));
end
OCMATFTE.Jpar=[OCMATFTE.dFDPAR;OCMATFTE.dFODPAR;OCMATFTE.dFEDPAR];
OCMATFTE.ODEcoord=1:size(OCMATFTE.Jpar,1);
 
OCMATFTE.Jext=[[OCMATFTE.dFDO OCMATFTE.dFDE];[OCMATFTE.dFODO OCMATFTE.dFODE];[OCMATFTE.dFEDO OCMATFTE.dFEDE]];

[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=1;

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