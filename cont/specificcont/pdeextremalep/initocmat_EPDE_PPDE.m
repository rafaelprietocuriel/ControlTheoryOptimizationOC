function [sol,opt]=initocmat_EPDE_PPDE(ppdeObj,ppdePrim,targetdistribution,opt,varargin)
% INITOCMAT_AE_EP initialization for asymptotic extremal calculation
%
% ppdeObj ...   corresponding optimal control model
% ppdePrim ...    equilibrium point
% contidx ... coordinates of the continuation variable
% targetdistribution ... determines direction of the continuation
%
% SOL=INITOCMAT_AE_EP(ppdeObj,ppdePrim,contidx,targetdistribution)


clear global OCMATCONT OCMATPPDESD
global OCMATCONT OCMATPPDESD
sol=[];
pathtype='';
truncationtime=[];
objectivevaluecalc=[];
if isempty(ppdeObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ppdePrim)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
if nargin==3
    opt=defaultocoptions;
end
if ~testconsistency(ppdePrim,ppdeObj,opt)
    ocmaterror('The elliptic PDE solution and oc model are not consistent.')
end
truncationtimeidx=find(strcmpi(varargin,'truncationtime'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(truncationtimeidx)
    truncationtime=varargin{truncationtimeidx+1};
end
if ~isempty(pathtypeidx)
    pathtype=varargin{pathtypeidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if isempty(pathtype)
    pathtype='s';
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
targetdistribution=targetdistribution(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)

if isempty(truncationtime)
    eval=eig(ppdePrim);
    eval(real(eval)>=0)=[];
    truncationtime=1/min(abs(real(eval)))*10;
end

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ppdeObj);
OCMATCONT.modelfunc=modelspecificfunc(ppdeObj,'4SaddlePathContinuation');

% initialize global variable (OCMATPPDESD) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATPPDESD.canonicalsystem=funch{1};
OCMATPPDESD.canonicalsystemjacobian=funch{2}{1};
OCMATPPDESD.canonicalsystemparameterjacobian=funch{2}{2};
OCMATPPDESD.canonicalsystemhessian=funch{3}{1};
OCMATPPDESD.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATPPDESD.bcinitial=funch{5}{1};
OCMATPPDESD.bcasymptotic=funch{5}{2};
OCMATPPDESD.bctransversality=funch{5}{3};

% function for Jacobian
OCMATPPDESD.bcjacobianinitial=funch{6}{1};
OCMATPPDESD.bcjacobianasymptotic=funch{6}{2};
OCMATPPDESD.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATPPDESD.geometry=funch{7}{1};
OCMATPPDESD.boundarycondition=funch{7}{2};
OCMATPPDESD.guard=funch{7}{3};
OCMATPPDESD.reset=funch{7}{4};
OCMATPPDESD.switchtime=funch{7}{5};
OCMATPPDESD.jacobianguard=funch{7}{7};
OCMATPPDESD.jacobianreset=funch{7}{8};
OCMATPPDESD.domaindiscretization=funch{7}{9};
OCMATPPDESD.timesettransformation=funch{7}{10};

% function for the calculation of the objective value during integration
OCMATPPDESD.objectivefunction=funch{8}{1};
OCMATPPDESD.objectivefunctionjacobian=funch{8}{2};
OCMATPPDESD.objectivefunctionparameterjacobian=funch{8}{3};
OCMATPPDESD.objectivefunctionderivativetime=funch{8}{4};

% general function
OCMATPPDESD.plotcontinuation=funch{11};
OCMATPPDESD.testadmissibility=funch{12};
OCMATPPDESD.datapath=funch{20};
OCMATPPDESD.saveintermediatefiles=funch{21};

% initialize solution using the repeated entries of the equilibrium
sol=generateodestruct(ppdePrim,TrivialArcMeshNum,truncationtime);

% add continuation parameter value
sol.parameters=[sol.parameters 0];

% reduce Jacobian to ODE part
J=jacobian(ppdePrim);
if issparse(J)
    J=full(J);
end
% mode and path specific variables
OCMATPPDESD.parametervalue=parametervalue(ppdeObj);
OCMATPPDESD.initialtime=sol.x0;
OCMATPPDESD.switchtimecoord=1:numel(sol.arcinterval)-2;
OCMATPPDESD.truncationtime=truncationtime;
OCMATPPDESD.saddlepoint=dependentvar(ppdePrim);
OCMATPPDESD.linearization=J;
OCMATPPDESD.pathtype=pathtype;
OCMATPPDESD.objectivevaluecalc=objectivevaluecalc;

ZeroDeviationTolerance=1e-5;
OCMATPPDESD.asymptoticmatrix=asymptoticbc(J,pathtype,'c',ZeroDeviationTolerance);
contcoordinate=1:numel(targetdistribution);
OCMATPPDESD.statecoordinate=contcoordinate;
OCMATPPDESD.costatecoordinate=contcoordinate(end)+contcoordinate;
OCMATPPDESD.totalcoordinate=1:size(sol.y,1);
OCMATPPDESD.statecostatenumcoordinate=1:2*statenum(ppdeObj);

OCMATPPDESD.targetcoordinate=contcoordinate;
OCMATPPDESD.startdistribution=OCMATPPDESD.saddlepoint(contcoordinate);
OCMATPPDESD.continuationvector=targetdistribution-OCMATPPDESD.startdistribution; % continue solution along the line from starting to target point
OCMATPPDESD.targetdistribution=targetdistribution;

%OCMATPPDESD.implicitcontrolindex=domaindata(statDistarcindex).implicitcontrolindex;
pathname=OCMATPPDESD.datapath();
[resultfile,globalvarfile]=OCMATPPDESD.saveintermediatefiles();
OCMATPPDESD.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATPPDESD.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATPPDESD.asymptoticmatrix,2);

OCMATCONT.codimension=1;

% space discretization
[OCMATPPDESD.geo OCMATPPDESD.bt]=csggeometry(ppdePrim);
OCMATPPDESD.points=points(ppdePrim);
OCMATPPDESD.edges=edges(ppdePrim);
OCMATPPDESD.triangles=triangles(ppdePrim);
OCMATPPDESD.boundarycondition=boundarycondition(ppdePrim);
[M invM bcG K Kadv]=femoperator(ppdePrim);
OCMATPPDESD.femop.M=sparse(M);
OCMATPPDESD.femop.invM=invM;
OCMATPPDESD.femop.bcG=bcG;
OCMATPPDESD.femop.K=sparse(K);
OCMATPPDESD.femop.Kadv=sparse(Kadv);
OCMATPPDESD.femop.invMKmKadv=sparse(M\(K-Kadv));
thresh=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
OCMATPPDESD.femop.inexactinvMKmKadv=sparse(OCMATPPDESD.femop.invMKmKadv.*(abs(OCMATPPDESD.femop.invMKmKadv)>thresh));
OCMATPPDESD.numpoints=size(OCMATPPDESD.points,2);
idx=1:length(sol.y(:,1));
N=length(idx)/OCMATPPDESD.numpoints;
OCMATPPDESD.coeffidx=reshape(idx,OCMATPPDESD.numpoints,N).';
OCMATPPDESD.movinghorizon=false;

if strcmp(getocoptions(opt,'GENERAL','BVPMethod'),'mtom0')
    MM=OCMATPPDESD.femop.M;
    MM(end+1:end+length(sol.parameters),end+1:end+length(sol.parameters))=eye(length(sol.parameters));
    opt=setocoptions(opt,'BVP','MassMatrix',MM,'LUsw',0);
    OCMATPPDESD.MassMatrix=true;
else
    OCMATPPDESD.MassMatrix=false;
end
np=size(OCMATPPDESD.points,2);
nt=size(OCMATPPDESD.triangles,2);

A=sparse(ones(3,1)*(1:nt),OCMATPPDESD.triangles(1:3,:),1,nt,np);
B=sparse(1:nt,1:nt,1./sum(A.'),nt,nt);
OCMATPPDESD.JacInt=(B*A).';
OCMATPPDESD.ppdeprimitive=ppdePrim;

if objectivevaluecalc
    OCMATPPDESD.trianglearea=pdetrg(OCMATPPDESD.points,OCMATPPDESD.triangles);
    OCMATPPDESD.triangleareajac=repmat(OCMATPPDESD.trianglearea,2*statenum(ppdeObj),1);
    t=truncationtime*sol.x;
    o=zeros(1,length(t));
    for ii=1:length(t)
        depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,sol.y(:,ii));
        tmp=OCMATPPDESD.objectivefunction(t(ii),OCMATPPDESD.points,depvarint,OCMATPPDESD.parametervalue,sol.arcarg);
        o(ii)=sum(OCMATPPDESD.trianglearea.*tmp);
    end
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(t))];
    OCMATPPDESD.objectivevaluecoord=size(sol.y,1);
end
