function sol= initocmat_FIMP_Popt(ocObj,hocTrj,optparindex,contindex,varargin)
%
%
% initocmat_FIMP_Popt
%
% if timeindex is empty, the endtime is assumed as the continuation parameter.
clear global OCMATCONT IOCMATFTE
global OCMATCONT IOCMATFTE
sol=[];
targetvalue=[];
targettype='';
initialcoordinate=[];
findoptimaljumptime=[];
phi=[];

if isempty(ocObj)
    ocmatmsg('impulse oc model is empty.')
    return
end
if isempty(hocTrj)
    ocmatmsg('hybrid oc trajectory is empty.')
    return
end

if ischar(optparindex)
    optparindex=parameterindex(ocObj,optparindex);
end

if isempty(optparindex)
    ocmatmsg('Parametervalue unknown.\n')
    return
end
targettypeidx=find(strcmpi(varargin,'targettype'));
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
initialcoordinateidx=find(strcmpi(varargin,'initialcoordinate'));
findoptimaljumptimeidx=find(strcmpi(varargin,'findoptimaljumptime'));
phiidx=find(strcmpi(varargin,'phi'));
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(initialcoordinateidx)
    initialcoordinate=varargin{initialcoordinateidx+1};
end
if ~isempty(findoptimaljumptimeidx)
    findoptimaljumptime=varargin{findoptimaljumptimeidx+1};
end
if ~isempty(phiidx)
    phi=varargin{phiidx+1};
end
if isempty(initialcoordinate)
    initialcoordinate=1:statenum(ocObj);
end
if isempty(findoptimaljumptime)
    findoptimaljumptime=false;
end
if isempty(targettype)
    targettype='parameter';
end
targettype=lower(targettype);

if strncmp(targettype,'p',1)
    contindex=parameterindex(ocObj,contindex);
    if isempty(contindex)
        return
    end
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
IOCMATFTE.objectivefunction=funch{8}{1};
IOCMATFTE.objectivefunctionjacobian=funch{8}{2};
IOCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
IOCMATFTE.salvagevalue=funch{5}{6};
IOCMATFTE.impulseobjectivefunction=funch{5}{7};
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
% mode and path specific variables
IOCMATFTE.parametervalue=parametervalue(ocObj);
IOCMATFTE.autonomous=isautonomous(ocObj);

sol=generatesolstruct(hocTrj,solver(hocTrj));
sol.parameters=[depvar([scoord cscoord],1).' depvar([scoord cscoord],end).'];
arctime=arcinterval(hocTrj);
jumparg=jumpargument(hocTrj);

switchtimeidx=find(jumparg>0 | [-1 jumparg(2:end-1) -1]==0);
if strncmp(targettype,'t',1)
    switchtimeidx=[switchtimeidx contindex];
end
sol.parameters=[sol.parameters IOCMATFTE.parametervalue(optparindex)];
sol.parameters=[sol.parameters arctime(switchtimeidx)];

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
    OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end


switch targettype(1)
    case 'p'
        sol.parameters=[sol.parameters IOCMATFTE.parametervalue(contindex)];
end
IOCMATFTE.initialtime=sol.x0;
IOCMATFTE.arctime=arctime;
IOCMATFTE.initialstate=depvar(initialcoordinate,1);
IOCMATFTE.initialdepvarcoord=[scoord(:).' cscoord(:).'];
IOCMATFTE.enddepvarcoord=2*snum+[scoord(:).' cscoord(:).'];
IOCMATFTE.optparametercoord=4*snum+1;
IOCMATFTE.switchtimecoord=IOCMATFTE.optparametercoord+(1:length(switchtimeidx));
IOCMATFTE.contcoord=IOCMATFTE.switchtimecoord(end)+1;
IOCMATFTE.switchtimeidx=switchtimeidx;
IOCMATFTE.varyarcintervalidx=zeros(1,length(arctime));
IOCMATFTE.varyarcintervalidx(switchtimeidx)=IOCMATFTE.switchtimecoord; % position of the free switching times in relation to the arcs
IOCMATFTE.initialcoordinate=initialcoordinate;
IOCMATFTE.contindex=contindex;
IOCMATFTE.targetvalue=targetvalue;
IOCMATFTE.jumparg=jumpargument(hocTrj);
IOCMATFTE.findoptimaljumptime=findoptimaljumptime;
IOCMATFTE.targettype=targettype;
IOCMATFTE.optparindex=optparindex;
IOCMATFTE.dOdp_phi=phi(:).';

if size(sol.y,1)==2*statenum(ocObj)
    OT=discountedsalvagevalue(ocObj,hocTrj)+impulseobjective(ocObj,hocTrj);
    o=objectivefunction(ocObj,hocTrj,1);
    sol.y(2*statenum(ocObj)+1,:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,hocTrj,1)))];
end
IOCMATFTE.objectivevaluecoord=size(sol.y,1);

pathname=IOCMATFTE.datapath();
[resultfile,globalvarfile]=IOCMATFTE.saveintermediatefiles();
IOCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
IOCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);

OCMATCONT.codimension=2;
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