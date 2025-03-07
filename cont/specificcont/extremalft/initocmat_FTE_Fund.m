function sol=initocmat_FTE_Fund(ocObj,ocTrj,contcoordinate,targetvalue,varargin)
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
fixendstate=[];
fixinitstate=[];
initialhorizon=[];
maxhorizon=[];
objectivevaluecalc=[];
hitvalue=[];
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
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(maxhorizon)
    maxhorizon=inf;
end
targetvalue=targetvalue(:);

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
    OCMATFTE.salvagevalue=funch{5}{6};
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

% test if pure state constraints are defined
OCMATFTE.stateconstraint=0;%stateconstraint(ocObj);
OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.maxhorizon=maxhorizon;
numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
end
if optimalhorizon
    sol.parameters=[sol.parameters initialhorizon];
end
OCMATFTE.switchtimecoord=1:numswitchtimes;
sol.parameters=[sol.parameters 0];
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+domaindata(ii).odedim^2;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder)+numel(domaindata(ii).daeorder)^2;%number of equations
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

% mode and path specific variables
depvar=dependentvar(ocTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);

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
OCMATFTE.targetcoordinate=contcoordinate;
if isempty(fixinitstate)
    OCMATFTE.startvalue=depvar(contcoordinate,1);
    OCMATFTE.continitstate=true;
elseif isempty(fixendstateidx)
    OCMATFTE.startvalue=depvar(contcoordinate,end);
    OCMATFTE.continitstate=false;
end
OCMATFTE.continuationvector=targetvalue-OCMATFTE.startvalue;
OCMATFTE.objectivevaluecalc=objectivevaluecalc;
OCMATFTE.optimalhorizon=optimalhorizon;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitstatecoord=fixinitstate;
if ~isempty(fixendstate)
    OCMATFTE.endstate=depvar(fixendstate,end);
    OCMATFTE.initstate=[];
elseif ~isempty(fixinitstate)
    OCMATFTE.initstate=depvar(fixinitstate,1);
    OCMATFTE.endstate=[];
end
pathname=OCMATFTE.datapath();
if objectivevaluecalc
    OCMATFTE.objectivevaluecoord=size(sol.y,1);
else
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(ii).numode,:);
end
if ~isempty(hitvalue)
    OCMATFTE.targetdistance=norm(OCMATFTE.startvalue-hitvalue)/norm(OCMATFTE.startvalue-targetvalue);
else
     OCMATFTE.targetdistance=[];
end
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);

OCMATCONT.codimension=1;

% variables for fundamental solution
n=statenum(ocObj);
OCMATFTE.solcoord=1:2*n;
OCMATFTE.fundcoord=2*n+(1:4*n^2);
OCMATFTE.fundsolinit=eye(2*n);
OCMATFTE.fundsolinit=OCMATFTE.fundsolinit(:);
OCMATFTE.fundmatcoord=reshape(1:4*n^2,2*n,[]);
OCMATFTE.fundnumcoord=4*n^2;

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