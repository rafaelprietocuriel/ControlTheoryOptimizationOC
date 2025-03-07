function sol=initocmat_INFTE_T(ocObj,ocTrj,timeindex,targetvalue,initialcoordinate,varargin)
%
% initocmat_FTE_T initialization for asymptotic extremal calculation
%
% SOL=initocmat_FTE_T(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
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
initialstate=[];
initialendtime=[];
initialarcargument=[];
fixendstate=[];
opt=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
optionidx=find(strcmpi(varargin,'option'));
initialstateidx=find(strcmpi(varargin,'initialstate'));
initialendtimeidx=find(strcmpi(varargin,'initialendtime'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end

%%% generate initial octrajectory if not provided
if isempty(ocTrj)
    if ~isempty(initialstateidx)
        initialstate=varargin{initialstateidx+1};
    end
    if ~isempty(initialendtimeidx)
        initialendtime=varargin{initialendtimeidx+1};
    end
    if ~isempty(initialarcargumentidx)
        initialarcargument=varargin{initialarcargumentidx+1};
    end
    if ~isempty(optionidx)
        opt=varargin{optionidx+1};
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    if isempty(initialarcargument)
        initialarcargument=0;
    end
    if isempty(initialendtime)
        initialendtime=0;
    end
    n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocTrj.x=linspace(0,1,n);
    initocPt.y=initialstate(:);
    initocPt.x=0;
    initocPt.arcarg=initialarcargument;
    initocPt.arcinterval=[0 initialendtime];
    initocPt.arcposition=[1;1];
    initialcostate=inftransversalitycondition(ocObj,octrajectory(initocPt));
    ocTrj.y=[initialstate(:);initialcostate];
    ocTrj.y=ocTrj.y(:,ones(1,n));
    ocTrj.arcposition=[1;n];
    ocTrj.arcinterval=[0 initialendtime];
    ocTrj.arcarg=initialarcargument;
    ocTrj.x0=0;
    ocTrj.timehorizon=initialendtime;
    ocTrj.modelparameter=parametervalue(ocObj);
    ocTrj.modelname=modelname(ocObj);
    ocTrj=octrajectory(ocTrj);
end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

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
OCMATFTE.bcinftransversality=funch{5}{3};

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

arctimes=arcinterval(ocTrj);

% test if pure state constraints are defined
OCMATFTE.stateconstraint=stateconstraint(ocObj);
numswitchtimes=length(sol.arcinterval(2:end));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=sol.arcinterval(2:end-1);
end
sol.parameters=[sol.parameters arctimes(timeindex)];
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
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);

OCMATFTE.initialtime=sol.x0;
OCMATFTE.initialstate=depvar(initialcoordinate,1);
OCMATFTE.switchtimecoord=1:numswitchtimes;
OCMATFTE.initialcoordinate=initialcoordinate;
OCMATFTE.varyparameterindex=timeindex;
OCMATFTE.targetvalue=targetvalue;

OCMATFTE.objectivevaluecalc=false;
if ~OCMATFTE.objectivevaluecalc
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(ii).numode,:);
end
OCMATFTE.fixendstatecoord=fixendstate;
if ~isempty(fixendstate)
    OCMATFTE.endstate=depvar(fixendstate,end);
end
pathname=OCMATFTE.datapath();
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