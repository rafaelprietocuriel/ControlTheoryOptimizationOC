function sol=initocmat_FIMP_RemI(ocObj,hocTrj,varargin)
%
%
% INITOCMAT_FIMP_ADDI add new impulse at time timeindex
clear global OCMATCONT IOCMATFTE
global OCMATCONT IOCMATFTE
sol=[];
initialcoordinate=[];
findoptimaljumptime=[];
objectivevaluecalc=[];
if isempty(ocObj)
    ocmatmsg('impulse oc model is empty.')
    return
end
if isempty(hocTrj)
    ocmatmsg('hybrid oc trajectory is empty.')
    return
end
initialcoordinateidx=find(strcmpi(varargin,'initialcoordinate'));
findoptimaljumptimeidx=find(strcmpi(varargin,'findoptimaljumptime'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(initialcoordinateidx)
    initialcoordinate=varargin{initialcoordinateidx+1};
end
if ~isempty(findoptimaljumptimeidx)
    findoptimaljumptime=varargin{findoptimaljumptimeidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end

if isempty(initialcoordinate)
    initialcoordinate=1:statenum(ocObj);
end
if isempty(findoptimaljumptime)
    findoptimaljumptime=false;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4FiniteImpulsePathContinuation');

% initialize global variable (IOCMATFTE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

IOCMATFTE.canonicalsystem=funch{1};
IOCMATFTE.canonicalsystemjacobian=funch{2}{1};
IOCMATFTE.canonicalsystemparameterjacobian=funch{2}{2};
%IOCMATFTE.canonicalsystemhessian=funch{3}{1};
%IOCMATFTE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
IOCMATFTE.bcinitial=funch{5}{1};
IOCMATFTE.bctransversality=funch{5}{2};
IOCMATFTE.bcoptimalhorizon=funch{5}{3};
IOCMATFTE.bcevent=funch{5}{4};
IOCMATFTE.bcinteriorevent=funch{5}{5};

% function for Jacobian
IOCMATFTE.bcjacobianinitial=funch{6}{1};
IOCMATFTE.bcjacobiantransversality=funch{6}{2};

% function describing the hybrid structure of the problem
IOCMATFTE.hybridinfo=funch{7}{1};
IOCMATFTE.domain=funch{7}{2};
IOCMATFTE.guard=funch{7}{3};
IOCMATFTE.reset=funch{7}{4};
IOCMATFTE.switchtime=funch{7}{5};
IOCMATFTE.jacobianguard=funch{7}{7};
IOCMATFTE.jacobianreset=funch{7}{8};
IOCMATFTE.domaindiscretization=funch{7}{9};
if objectivevaluecalc
    IOCMATFTE.objectivefunction=funch{8}{1};
    IOCMATFTE.objectivefunctionjacobian=funch{8}{2};
    IOCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
end
if ~isautonomous(ocObj)
    IOCMATFTE.canonicalsystemderivativetime=funch{2}{3};
    IOCMATFTE.objectivefunctionderivativetime=funch{8}{4};
end

% general function
IOCMATFTE.plotcontinuation=funch{11};
IOCMATFTE.testadmissibility=funch{12};
IOCMATFTE.datapath=funch{20};
IOCMATFTE.saveintermediatefiles=funch{21};

hybridinfo=IOCMATFTE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=IOCMATFTE.domain(hybridinfo.arcarg(ii));
end
depvar=dependentvar(hocTrj);
snum=statenum(ocObj);
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);

sol=generatesolstruct(hocTrj,solver(hocTrj));
sol.parameters=[depvar(:,1).' depvar(:,end).'];
arctime=arcinterval(hocTrj);
jumparg=jumpargument(hocTrj);
fixedtimeidx=find(jumparg<0);
switchtimeidx=find(jumparg>0 | [-1 jumparg(2:end-1) -1]==0);

sol.parameters=[sol.parameters arctime(switchtimeidx) 1];
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
IOCMATFTE.parametervalue=parametervalue(ocObj);
IOCMATFTE.autonomous=isautonomous(ocObj);

IOCMATFTE.initialtime=sol.x0;
IOCMATFTE.endtime=arctime(end);
IOCMATFTE.arctime=arctime;
IOCMATFTE.initialstate=depvar(initialcoordinate,1);
IOCMATFTE.initialdepvarcoord=[scoord(:).' cscoord(:).'];
IOCMATFTE.enddepvarcoord=2*snum+[scoord(:).' cscoord(:).'];
IOCMATFTE.switchtimecoord=4*snum+[1:length(switchtimeidx)];
IOCMATFTE.fixedtimeidx=fixedtimeidx;
IOCMATFTE.switchtimeidx=switchtimeidx;
IOCMATFTE.initialcoordinate=initialcoordinate;
IOCMATFTE.jumparg=jumparg;
IOCMATFTE.targetvalue=0;
IOCMATFTE.findoptimaljumptime=findoptimaljumptime;
IOCMATFTE.objectivevaluecalc=objectivevaluecalc;

if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,hocTrj,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,hocTrj,1)))];
end
if objectivevaluecalc
    IOCMATFTE.objectivevaluecoord=size(sol.y,1);
else
    IOCMATFTE.objectivevaluecoord=[];
end

if ~IOCMATFTE.objectivevaluecalc
    sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(ii).numode,:);
    IOCMATFTE.objectivevaluecoord=[];
end
IOCMATFTE.varyendtime=false;

pathname=IOCMATFTE.datapath();
[resultfile,globalvarfile]=IOCMATFTE.saveintermediatefiles();
IOCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
IOCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);

OCMATCONT.codimension=1;
OCMATCONT.TargetValueNum=1;

function sol=generatesolstruct(hocTrj,solvername,varargin)

sol.x=independentvar(hocTrj);
sol.x([1 end])=[];
sol.y=dependentvar(hocTrj);
sol.y(:,[1 end])=[];
sol.arcarg=arcargument(hocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(hocTrj);
sol.arcposition=arcposition(hocTrj)-1;
sol.solver=solvername;

if isfield(hocTrj.solverinfo,'tangent')
    sol.solverinfo.tangent=hocTrj.solverinfo.tangent;
else
    sol.solverinfo.tangent=[];
end
if isfield(hocTrj.solverinfo,'coeff')
    sol.solverinfo.coeff=hocTrj.solverinfo.coeff;
else
    sol.solverinfo.coeff=[];
end
if isfield(hocTrj.solverinfo,'yp')
    sol.solverinfo.yp=hocTrj.solverinfo.yp;
end
if isfield(hocTrj.solverinfo,'ypmid')
    sol.solverinfo.ypmid=hocTrj.solverinfo.ypmid;
end
if isfield(hocTrj.solverinfo,'ymid')
    sol.solverinfo.yp=hocTrj.solverinfo.yp;
end
if isfield(hocTrj.solverinfo,'ymid')
    sol.solverinfo.ymid=hocTrj.solverinfo.ymid;
end