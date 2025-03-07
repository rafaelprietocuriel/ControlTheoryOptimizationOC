function sol=initocmat_USER_AE(ocObj,ocTrj,contcoordinate,targetvalue,varargin)
%

clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
freehorizon=[];
objectivevaluecalc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('oc trajectory is empty.')
    return
end
freehorizonidx=find(strcmpi(varargin,'freehorizon'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(freehorizonidx)
    freehorizon=varargin{freehorizonidx+1};
end
if isempty(freehorizon)
    freehorizon=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end

targetvalue=targetvalue(:);

% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4UserContinuation');

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
OCMATAE.bcend=funch{5}{2};
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
if objectivevaluecalc
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
end

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};
OCMATAE.initializesol=funch{22};
hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

sol=OCMATAE.initializesol(ocTrj);
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
numode=OCMATCONT.DOMAINDDATA(1).numeq;
% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.parameters);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;

OCMATAE.objectivevaluecalc=objectivevaluecalc;
OCMATAE.freehorizon=freehorizon;
if freehorizon
    if isempty(OCMATAE.switchtimecoord)
        OCMATAE.freehorizoncoord=1;
    else
        OCMATAE.freehorizoncoord=length(sol.parameters)+1;
    end
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
sol.parameters=[sol.parameters 0];
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;

if objectivevaluecalc && length(sol.y(:,1))==numode-1
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocTrj,1)))];
end
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=numode;
end
