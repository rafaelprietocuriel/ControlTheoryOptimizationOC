function [sol opt]=initocmat_PPDE_PPDE(ppdeObj,ppdeAsym,targetdistribution,opt,varargin)
%
% INITOCMAT_AE_AE initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_AE(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_AE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATPPDESD
global OCMATCONT OCMATPPDESD
sol=[];
pathtpe=pathtype(ppdeAsym);
objectivevaluecalc=[];
movinghorizon=[];
targetvalue=[];
targettype='';
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
if isempty(pathtpe)
    pathtpe='s';
end
if isempty(ppdeObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ppdeAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==3
    opt=defaultocoptions;
end
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end

if isfield(ppdeAsym.discretizationinfo.coord,'objectivevalue') && ~isempty(ppdeAsym.discretizationinfo.coord.objectivevalue)
    oldobjectivevaluecalc=1;
    objectivevaluecoord=ppdeAsym.discretizationinfo.coord.objectivevalue;
else
    oldobjectivevaluecalc=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=oldobjectivevaluecalc;
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if isempty(movinghorizon)
    movinghorizon=0;
end

targetdistribution=targetdistribution(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ppdeObj);
OCMATCONT.modelfunc=modelspecificfunc(ppdeObj,'4SaddlePathContinuation');

% initialize global variable (OCMATPPDESD) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
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

sol=generatesolstruct(ppdeAsym.ppdetrajectory,solver(ppdeAsym));

% test if pure state constraints are defined
sol.parameters=[sol.parameters 0];
% mode and path specific variables
ppdePrim=ppdeAsym.limitset;
J=linearization(ppdePrim);
depvar=dependentvar(ppdeAsym);
OCMATPPDESD.parametervalue=parametervalue(ppdeObj);

OCMATPPDESD.initialtime=ppdeAsym.discretizationinfo.timeinterval(1);
switch pathtpe
    case {'s','sc','cs'}
        OCMATPPDESD.truncationtime=ppdeAsym.discretizationinfo.timeinterval(end);
    case 'u'
        OCMATPPDESD.truncationtime=-ppdeAsym.discretizationinfo.timeinterval(end);
end
OCMATPPDESD.linearization=J;
OCMATPPDESD.saddlepoint=dependentvar(ppdePrim);
ZeroDeviationTolerance=1e-5;
OCMATPPDESD.asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
contcoordinate=1:numel(targetdistribution);

OCMATPPDESD.targetcoordinate=contcoordinate;
OCMATPPDESD.statecoordinate=contcoordinate;
OCMATPPDESD.startdistribution=depvar(contcoordinate,1);
OCMATPPDESD.continuationvector=targetdistribution-OCMATPPDESD.startdistribution;
OCMATPPDESD.totalcoordinate=ppdeAsym.discretizationinfo.coord.totaldepvar;
OCMATPPDESD.targetdistribution=targetdistribution;

OCMATPPDESD.pathtype=pathtpe;

pathname=OCMATPPDESD.datapath();
[resultfile,globalvarfile]=OCMATPPDESD.saveintermediatefiles();
OCMATPPDESD.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATPPDESD.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATPPDESD.switchtimecoord=[];
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
OCMATPPDESD.coeffidx=ppdeAsym.discretizationinfo.coord.meshidx;
OCMATPPDESD.ppdeprimitive=ppdePrim;

np=size(OCMATPPDESD.points,2);
nt=size(OCMATPPDESD.triangles,2);

A=sparse(ones(3,1)*(1:nt),OCMATPPDESD.triangles(1:3,:),1,nt,np);
B=sparse(1:nt,1:nt,1./sum(A.'),nt,nt);
OCMATPPDESD.JacInt=(B*A).';

if objectivevaluecalc
    OCMATPPDESD.trianglearea=pdetrg(OCMATPPDESD.points,OCMATPPDESD.triangles);
    OCMATPPDESD.triangleareajac=repmat(OCMATPPDESD.trianglearea,2*statenum(ppdeObj),1);
    if ~oldobjectivevaluecalc
        t=truncationtime*sol.x;
        o=zeros(1,length(t));
        for ii=1:length(t)
            depvarint=pdeintrp(OCMATPPDESD.edges,OCMATPPDESD.triangles,sol.y(:,ii));
            tmp=truncationtime*OCMATPPDESD.objectivefunction(t(ii),OCMATPPDESD.points,depvarint,OCMATPPDESD.parametervalue,sol.arcarg);
            o(ii)=sum(OCMATPPDESD.trianglearea.*tmp);
        end
        sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(t))];
    end
    OCMATPPDESD.objectivevaluecoord=size(sol.y,1);
elseif oldobjectivevaluecalc
    sol.y(oldobjectivevaluecoord,:)=[];
end
OCMATPPDESD.movinghorizon=movinghorizon;
if movinghorizon
    OCMATPPDESD.distance=norm(OCMATPPDESD.saddlepoint-depvar(OCMATPPDESD.totalcoordinate,end));
    sol.parameters=[sol.timeinterval(end) sol.parameters];
    OCMATPPDESD.movinghorizoncoord=1;
end

if strcmp(getocoptions(opt,'GENERAL','BVPMethod'),'mtom0')
    MM=OCMATPPDESD.femop.M;
    MM(end+1:end+length(sol.parameters),end+1:end+length(sol.parameters))=eye(length(sol.parameters));
    opt=setocoptions(opt,'BVP','MassMatrix',MM);
    OCMATPPDESD.MassMatrix=true;
else
    OCMATPPDESD.MassMatrix=false;
end

OCMATPPDESD.objectivevaluecalc=objectivevaluecalc;

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.timeinterval=timeinterval(ocTrj);
sol.arcarg=0;
sol.parameters=freeparameters(ocTrj);
sol.solver=solvername;

if isfield(ocTrj.discretizationinfo,'tangent')
    sol.discretizationinfo.tangent=ocTrj.discretizationinfo.tangent;
else
    sol.discretizationinfo.tangent=[];
end
if isfield(ocTrj.discretizationinfo,'coeff')
    sol.discretizationinfo.coeff=ocTrj.discretizationinfo.coeff;
else
    sol.discretizationinfo.coeff=[];
end
if isfield(ocTrj.discretizationinfo,'yp')
    sol.discretizationinfo.yp=ocTrj.discretizationinfo.yp;
end
if isfield(ocTrj.discretizationinfo,'ypmid')
    sol.discretizationinfo.ypmid=ocTrj.discretizationinfo.ypmid;
end
if isfield(ocTrj.discretizationinfo,'ymid')
    sol.discretizationinfo.yp=ocTrj.discretizationinfo.yp;
end
if isfield(ocTrj.discretizationinfo,'ymid')
    sol.discretizationinfo.ymid=ocTrj.discretizationinfo.ymid;
end