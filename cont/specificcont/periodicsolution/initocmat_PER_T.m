function sol=initocmat_PER_T(ocObj,ocPS,timeindex,varargin)
% INITOCMAT_PS_P initialization for the continuation of a peridoic solution
% with respect to a parameter
%
% SOL=INITOCMAT_PS_P(OCOBJ,OCPS,PARINDEX)
%
% SOL=INITOCMAT_PS_P(OCOBJ,OCPS,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATPS
global OCMATCONT OCMATPS
sol=[];
targettime=[];
linearizationcalc=[];
varyperiod=[];
objectivevaluecalc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocPS)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin>=4
    targettime=varargin{1};
end
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
linearizationidx=find(strcmpi(varargin,'linearization'));
varyperiodidx=find(strcmpi(varargin,'varyperiod'));
if ~isempty(linearizationidx)
    linearizationcalc=varargin{linearizationidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if isempty(linearizationcalc)
    linearizationcalc=0;
end
if ~isempty(varyperiodidx)
    varyperiod=varargin{varyperiodidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4PeriodicSolutionContinuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATPS.canonicalsystem=funch{1};
OCMATPS.canonicalsystemjacobian=funch{2}{1};
OCMATPS.canonicalsystemparameterjacobian=funch{2}{2};
OCMATPS.canonicalsystemhessian=funch{3}{1};
OCMATPS.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATPS.bcperiodic=funch{5}{1};

% function for Jacobian
OCMATPS.bcjacobianperiodic=funch{6}{1};

% function describing the hybrid structure of the problem
OCMATPS.hybridinfo=funch{7}{1};
OCMATPS.domain=funch{7}{2};
OCMATPS.guard=funch{7}{3};
OCMATPS.reset=funch{7}{4};
OCMATPS.switchtime=funch{7}{5};
OCMATPS.jacobianguard=funch{7}{7};
OCMATPS.jacobianreset=funch{7}{8};
OCMATPS.domaindiscretization=funch{7}{9};
OCMATPS.timesettransformation=funch{7}{10};
if ~isautonomous(ocObj)
    OCMATPS.canonicalsystemderivativetime=funch{2}{3};
    OCMATPS.objectivefunction=funch{8}{1};
    OCMATPS.objectivefunctionjacobian=funch{8}{2};
    OCMATPS.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATPS.objectivefunctionderivativetime=funch{8}{4};
end

% general function
OCMATPS.plotcontinuation=funch{11};
OCMATPS.testadmissibility=funch{12};
OCMATPS.datapath=funch{20};
OCMATPS.saveintermediatefiles=funch{21};

hybridinfo=OCMATPS.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATPS.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocPS,solver(ocPS));
if isfield(sol,'parameters')
    sol.parameters=[sol.parameters arcinterval(ocPS,timeindex)];
else
    sol.parameters=arcinterval(ocPS,timeindex);
end
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivevalue(ocObj,ocPS,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocPS,1)))];
    OCMATPS.objectivevaluecoord=size(sol.y,1);
end
if ~objectivevaluecalc && length(sol.y(:,1))>domaindata(ii).odedim
    sol.y=sol.y(1:numode,:);
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
OCMATPS.varyperiod=varyperiod;
if OCMATPS.varyperiod
    OCMATPS.periodfunc=funch{9};
else
    OCMATPS.periodfunc=[];
end
    
% mode and path specific variables
OCMATPS.parametervalue=parametervalue(ocObj);
OCMATPS.autonomous=isautonomous(ocObj);
OCMATPS.initialtime=sol.x0;
OCMATPS.period=period(ocPS);
OCMATPS.switchtimecoord=1:numel(sol.arcinterval)-2;
OCMATPS.varytimeindex=timeindex;
OCMATPS.targettime=targettime;
OCMATPS.linearizationcalc=linearizationcalc;
OCMATPS.objectivevaluecalc=objectivevaluecalc;
if linearizationcalc
    OCMATPS.linearizationinit=eye(2*statenum(ocObj));
    OCMATPS.linearizationinit=OCMATPS.linearizationinit(:);
    OCMATPS.linearizationindex=2*statenum(ocObj)+[1:4*statenum(ocObj)^2];
    OCMATPS.linearizationindex=reshape(OCMATPS.linearizationindex,2*statenum(ocObj),2*statenum(ocObj));
end
OCMATCONT.codimension=1;

pathname=OCMATPS.datapath();
[resultfile,globalvarfile]=OCMATPS.saveintermediatefiles();
OCMATPS.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATPS.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATCONT.numintegralconstraint=0;
