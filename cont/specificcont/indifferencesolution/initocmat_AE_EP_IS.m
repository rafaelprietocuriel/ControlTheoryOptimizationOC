function sol=initocmat_AE_EP_IS(ocObj,ocAsym,parindex,opt,varargin)
% initocmat_AE_EP_IS initialization for an indifferencesolution at a
% (stalling) equilibrium


clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];
pathtpe=pathtype(ocAsym);
targetparametervalue=[];
targetparameterindex=[];
freeendtime=[];
objectivevaluecalc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
stableflag=~isempty(strfind(pathtpe,'s'));
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(parindex)
    ocmatmsg('No or unknown parameter specified.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
targetparameterindexidx=find(strcmpi(varargin,'targetparameterindex'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(targetparameterindexidx)
    targetparameterindex=varargin{targetparameterindexidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if isempty(freeendtime)
    freeendtime=false;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(targetparameterindex) && ~isempty(targetparametervalue)
    targetparameterindex=parindex(1);
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATINDIF.canonicalsystem=funch{1};
OCMATINDIF.canonicalsystemjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
OCMATINDIF.canonicalsystemhessian=funch{3}{1};
OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.equilibrium=funch{5}{4};
OCMATINDIF.indifferenceequilibrium=funch{5}{5};
% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcjacobianasymptotic=funch{6}{2};
OCMATINDIF.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATINDIF.hybridinfo=funch{7}{1};
OCMATINDIF.domain=funch{7}{2};
OCMATINDIF.guard=funch{7}{3};
OCMATINDIF.reset=funch{7}{4};
OCMATINDIF.switchtime=funch{7}{5};
OCMATINDIF.jacobianguard=funch{7}{7};
OCMATINDIF.jacobianreset=funch{7}{8};
OCMATINDIF.domaindiscretization=funch{7}{9};
OCMATINDIF.timesettransformation=funch{7}{10};
if objectivevaluecalc
    OCMATINDIF.objectivefunction=funch{8}{1};
    OCMATINDIF.objectivefunctionjacobian=funch{8}{2};
    OCMATINDIF.objectivefunctionparameterjacobian=funch{8}{3};
end
% general function
OCMATINDIF.plotcontinuation=funch{11};
OCMATINDIF.testadmissibility=funch{12};
OCMATINDIF.datapath=funch{20};
OCMATINDIF.saveintermediatefiles=funch{21};

hybridinfo=OCMATINDIF.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATINDIF.domain(hybridinfo.arcarg(ii));
end

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(1).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% mode and path specific variables
limSet=limitset(ocAsym);
limSetdepvar=dependentvar(limSet);
limSetDim=length(limSetdepvar);

depvar=dependentvar(ocAsym);
JacobianMatrix=linearization(limSet);

OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=initialtime(ocAsym);
%OCMATINDIF.switchtimecoord=1:numel(sol.arcinterval)-2;
OCMATINDIF.inftimetransformation=inftimetransformation(ocAsym);
OCMATINDIF.linearization=linearization(limSet);
OCMATINDIF.saddlepoint=dependentvar(limSet);
OCMATINDIF.pathtype=pathtpe;
OCMATINDIF.objectivevaluecalc=objectivevaluecalc;
OCMATINDIF.stableflag=stableflag;
[dum OCMATINDIF.numstable OCMATINDIF.numunstable OCMATINDIF.numcenter]=asymptoticbc(JacobianMatrix,pathtpe,'c',ZeroDeviationTolerance);
OCMATINDIF.freeparameterindex=parindex;
OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.targetparameterindex=targetparameterindex;

OCMATINDIF.freeendtime=freeendtime;
OCMATCONT.codimension=1;

pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.equilibriumcoord=1:length(limSetdepvar);
switch pathtpe
    case 's'
        OCMATINDIF.subspacedim=OCMATINDIF.numstable;
    case 'u'
        OCMATINDIF.subspacedim=OCMATINDIF.numunstable;
    case {'sc','cs'}
        OCMATINDIF.subspacedim=OCMATINDIF.numstable+OCMATINDIF.numcenter;
end
OCMATINDIF.orthspacedim=limSetDim-OCMATINDIF.subspacedim;
Y=zeros(OCMATINDIF.orthspacedim,OCMATINDIF.subspacedim);
OCMATINDIF.Y=Y;
% orthogonal basis for stable eigenspace
OCMATINDIF.Q0=computeBase(JacobianMatrix,stableflag,OCMATINDIF.subspacedim);
OCMATCONT.HE.Ycoord=reshape(length(limSetdepvar)+(1:OCMATINDIF.orthspacedim*OCMATINDIF.subspacedim),OCMATINDIF.orthspacedim,OCMATINDIF.subspacedim);
OCMATINDIF.Id=eye(OCMATINDIF.orthspacedim);
OCMATINDIF.numY=numel(Y);

sol=generatesolstruct(ocAsym);
numswitchtimes=length(sol.arcinterval(2:end-1));
sol.parameters=[sol.arcinterval(2:end-1)];
sol.parameters=[sol.parameters limSetdepvar.' Y(:).'];
if freeendtime
    yend=depvar(:,end);
    hatx=OCMATINDIF.saddlepoint;
    OCMATINDIF.distance=norm(yend-hatx);
    OCMATINDIF.freeendtimecoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
OCMATINDIF.freeparametercoordinate=length(sol.parameters)+(1:length(parindex));
sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];
OCMATINDIF.switchtimecoord=1:numswitchtimes;
OCMATINDIF.truncationtime=sol.arcinterval(end);
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocAsym,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocAsym,1)))];
end
if objectivevaluecalc
    OCMATINDIF.objectivevaluecoord=size(sol.y,1);
end

OCMATCONT.HE.numinitialcondition=1;
OCMATCONT.HE.numendcondition=OCMATINDIF.orthspacedim;
OCMATCONT.HE.Ycoord=OCMATCONT.HE.Ycoord+numel(sol.arcarg)-1;
OCMATCONT.HE.equilibriumcoord=OCMATCONT.HE.equilibriumcoord+numel(sol.arcarg)-1;

function  Q=computeBase(J,stableflag,NSub)

if ~stableflag
    J=-J;
end
%[VU, DU] = realeig(J);
% Select first NSub eigenvectors: unstable eigenspace
% Compute orthonormal basis for the eigenspace
%VU = VU(:,1:NSub);
%[Q,RU] = qr(VU);

[Qt,RUt]=schur(J);
[Q,TS]=ordschur(Qt,RUt,'lhp');
%[Q,S,V]=svd(Q(:,1:NSub));




function sol=generatesolstruct(ocTrj,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
sol.solver='';%solvername;

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

