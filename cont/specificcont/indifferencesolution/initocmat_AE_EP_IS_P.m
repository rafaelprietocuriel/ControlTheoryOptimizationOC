function sol=initocmat_AE_EP_IS_P(ocObj,ocAsym,parindex,initialcoordinate,opt,varargin)
% INITOCMAT_AEP_EP initialization for asymptotic extremal calculation
% varying a parameter
%
% SOL=INITOCMAT_AEP_EP(ocObj,ocEP,contidx,targetvalue)


clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtpe=pathtype(ocAsym);
targetparametervalue=[];
targetfunction=[];
followequilibrium=[];
movinghorizon=[];
objectivevaluecalc=[];
findoptimalparameter=[];
jacobian=[];
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
fixcoordinateidx=find(strcmpi(varargin,'fixcoordinate'));
jacobianidx=find(strcmpi(varargin,'jacobian'));
arcid4epidx=find(strcmpi(varargin,'arcid4ep'));
excludecoordinate4epidx=find(strcmpi(varargin,'excludecoordinate4ep'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
targetfunctionidx=find(strcmpi(varargin,'targetfunction'));
followequilibriumidx=find(strcmpi(varargin,'followequilibrium'));
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
findoptimalparameteridx=find(strcmpi(varargin,'findoptimalparameter'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(fixcoordinateidx)
    fixcoordinate=varargin{fixcoordinateidx+1};
else
    fixcoordinate=[];
end
if ~isempty(jacobianidx)
    jacobian=varargin{jacobianidx+1};
end
if ~isempty(excludecoordinate4epidx)
    excludecoordinate4ep=varargin{excludecoordinate4epidx+1};
else
    excludecoordinate4ep=[];
end
if ~isempty(arcid4epidx)
    arcid4ep=varargin{arcid4epidx+1};
else
    arcid4ep=[];
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(targetfunctionidx)
    targetfunction=varargin{targetfunctionidx+1};
end
if ~isempty(followequilibriumidx)
    followequilibrium=varargin{followequilibriumidx+1};
end
if ~isempty(findoptimalparameteridx)
    findoptimalparameter=varargin{findoptimalparameteridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if isempty(movinghorizon)
    movinghorizon=false;
end
if isempty(followequilibrium)
    followequilibrium=false;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(findoptimalparameter)
    findoptimalparameter=false;
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.equilibrium=funch{5}{4};

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
if findoptimalparameter
    OCMATAE.hamiltonianderivative=funch{9};
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};
if ~isempty(targetfunction)
    OCMATAE.targetfunction=funch{22};
else
    OCMATAE.targetfunction=[];
end
hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
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
limSetdepvar(excludecoordinate4ep)=[];
limSetDim=length(limSetdepvar);

depvar=dependentvar(ocAsym);
if ~isempty(jacobian)
    JacobianMatrix=jacobian;
else
    JacobianMatrix=linearization(limSet);
    JacobianMatrix(:,excludecoordinate4ep)=[];
    JacobianMatrix(excludecoordinate4ep,:)=[];
end
if ~isempty(arcid4ep)
    OCMATAE.arcid4ep=arcid4ep;
else
    OCMATAE.arcid4ep=arcargument(limSet);
end
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=initialtime(ocAsym);
%OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;
OCMATAE.inftimetransformation=inftimetransformation(ocAsym);
OCMATAE.linearization=linearization(limSet);
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.saddlepoint(excludecoordinate4ep)=[];
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=objectivevaluecalc;
OCMATAE.stableflag=stableflag;
[dum OCMATAE.numstable OCMATAE.numunstable OCMATAE.numcenter]=asymptoticbc(JacobianMatrix,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.initialstate=depvar(initialcoordinate,1);
OCMATAE.initialcoordinate=initialcoordinate;
OCMATAE.freeparameterindex=parindex;
OCMATAE.targetparametervalue=targetparametervalue;
OCMATAE.fixcoordinate=fixcoordinate;
OCMATAE.excludecoordinate4ep=excludecoordinate4ep;
OCMATAE.varcoordinate=1:length(limSetdepvar);
OCMATAE.followequilibrium=followequilibrium;
if ~isempty(fixcoordinate)
    OCMATAE.fixcoordinatevalue=limSetdepvar(OCMATAE.fixcoordinate);
    limSetdepvar(OCMATAE.fixcoordinate)=[];
    OCMATAE.varcoordinate(OCMATAE.fixcoordinate)=[];
end
OCMATAE.movinghorizon=movinghorizon;
OCMATCONT.codimension=1;
OCMATAE.findoptimalparameter=findoptimalparameter;

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.equilibriumcoord=1:length(limSetdepvar);
switch pathtpe
    case 's'
        OCMATAE.subspacedim=OCMATAE.numstable;
    case 'u'
        OCMATAE.subspacedim=OCMATAE.numunstable;
    case {'sc','cs'}
        OCMATAE.subspacedim=OCMATAE.numstable+OCMATAE.numcenter;
end
OCMATAE.orthspacedim=limSetDim-OCMATAE.subspacedim;
Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);
OCMATAE.Y=Y;
% orthogonal basis for stable eigenspace
OCMATAE.Q0=computeBase(JacobianMatrix,stableflag,OCMATAE.subspacedim);
OCMATCONT.HE.Ycoord=reshape(length(limSetdepvar)+(1:OCMATAE.orthspacedim*OCMATAE.subspacedim),OCMATAE.orthspacedim,OCMATAE.subspacedim);
OCMATAE.Id=eye(OCMATAE.orthspacedim);
OCMATAE.numY=numel(Y);

sol=generatesolstruct(ocAsym);
numswitchtimes=length(sol.arcinterval(2:end-1));
sol.parameters=[sol.arcinterval(2:end-1)];
sol.parameters=[sol.parameters limSetdepvar.' Y(:).'];
if movinghorizon
    yend=depvar(:,end);
    hatx=OCMATAE.saddlepoint;
    if ~isempty(OCMATAE.excludecoordinate4ep)
        yend(OCMATAE.excludecoordinate4ep,:)=[];
    end
    OCMATAE.distance=norm(yend-hatx);
    OCMATAE.movinghorizoncoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
OCMATAE.freeparametercoordinate=length(sol.parameters)+(1:length(parindex));
sol.parameters=[sol.parameters OCMATAE.parametervalue(parindex)];
OCMATAE.switchtimecoord=1:numswitchtimes;
OCMATAE.truncationtime=sol.arcinterval(end);
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocAsym,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocAsym,1)))];
end
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=size(sol.y,1);
end

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);
OCMATCONT.HE.numendcondition=OCMATAE.orthspacedim;
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

