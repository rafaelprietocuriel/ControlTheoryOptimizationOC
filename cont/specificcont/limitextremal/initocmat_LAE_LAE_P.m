function sol=initocmat_LAE_LAE_P(ocObj,ocAsym,parindex,opt,varargin)
% INITOCMAT_LAE_LAE initialization for limit point solution continuation
% contInfo is the structure returned at the limit point solution

clear global OCMATCONT OCMATLSC
global OCMATCONT OCMATLSC
sol=[];
targetvalue=[];
freevector=[];
targetcoordinate=[];
targetparametervalue=[];
objectivevaluecalc=[];
freeendtime=[];
pathtpe=pathtype(ocAsym);
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
if nargin==3 
    opt=defaultocoptions;
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
freevectoridx=find(strcmpi(varargin,'freevector'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if isempty(freeendtime)
    freeendtime=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
limSet=limitset(ocAsym);
if isperiodic(limSet)
    limitsettype='d';
    if strcmp(pathtpe,'s')
        freeendtime=-1;
    end
else
    limitsettype='c';
end

OCMATCONT.OPTIONS.SymDerivative=2;
BVPMethod=getocoptions(opt,'GENERAL','BVPMethod');
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

OCMATCONT.codimension=2;
OCMATLSC.stateconstraint=0;

%OCMATLSC.freeendtime=[];
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
%OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4LimitPathContinuation');
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATLSC) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATLSC.canonicalsystem=funch{1};
OCMATLSC.canonicalsystemjacobian=funch{2}{1};
OCMATLSC.canonicalsystemparameterjacobian=funch{2}{2};
OCMATLSC.canonicalsystemhessian=funch{3}{1};
OCMATLSC.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATLSC.bcinitial=funch{5}{1};
OCMATLSC.bcasymptotic=funch{5}{2};
OCMATLSC.bctransversality=funch{5}{3};
OCMATLSC.equilibrium=funch{5}{4};

% function for Jacobian
OCMATLSC.bcjacobianinitial=funch{6}{1};
OCMATLSC.bcjacobianasymptotic=funch{6}{2};
OCMATLSC.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATLSC.hybridinfo=funch{7}{1};
OCMATLSC.domain=funch{7}{2};
OCMATLSC.guard=funch{7}{3};
OCMATLSC.reset=funch{7}{4};
OCMATLSC.switchtime=funch{7}{5};
OCMATLSC.jacobianguard=funch{7}{7};
OCMATLSC.jacobianreset=funch{7}{8};
OCMATLSC.domaindiscretization=funch{7}{9};
OCMATLSC.timesettransformation=funch{7}{10};

% general function
OCMATLSC.plotcontinuation=funch{11};
OCMATLSC.testadmissibility=funch{12};
OCMATLSC.datapath=funch{20};
OCMATLSC.saveintermediatefiles=funch{21};

hybridinfo=OCMATLSC.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATLSC.domain(hybridinfo.arcarg(ii));
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
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATLSC.statecoord=scoord(:).';
OCMATLSC.statecostatecoord=[scoord(:).' cscoord(:).'];

OCMATLSC.freevector=freevector;
OCMATLSC.targetvalue=targetvalue;
OCMATLSC.targetcoordinate=targetcoordinate;
OCMATLSC.varyparameterindex=parindex;
OCMATLSC.targetparametervalue=targetparametervalue;

depvar=dependentvar(ocAsym);

sol=generatesolstruct(ocAsym,BVPMethod);

sol.parameters=sol.arcinterval(2:end-1);
OCMATLSC.switchtimecoord=1:length(sol.parameters);


% mode and path specific variables
limSet=limitset(ocAsym);
limSetdepvar=dependentvar(limSet);
limSetDim=length(limSetdepvar);
JacobianMatrix=linearization(limSet);
OCMATLSC.parametervalue=parametervalue(ocObj);
OCMATLSC.initialtime=sol.x0;
OCMATLSC.startvalue=depvar(scoord,1);

OCMATLSC.inftimetransformation=inftimetransformation(ocAsym);
OCMATLSC.truncationtime=sol.arcinterval(end);
OCMATLSC.linearization=JacobianMatrix;
OCMATLSC.saddlepoint=dependentvar(limSet);
OCMATLSC.saddlepoint=OCMATLSC.saddlepoint(:,1);
[OCMATLSC.asymptoticmatrix OCMATLSC.numstable OCMATLSC.numunstable OCMATLSC.numcenter]=asymptoticbc(JacobianMatrix,pathtpe,limitsettype,ZeroDeviationTolerance);
OCMATLSC.pathtype=pathtpe;
OCMATLSC.stableflag=stableflag;

pathname=OCMATLSC.datapath();
[resultfile,globalvarfile]=OCMATLSC.saveintermediatefiles();
OCMATLSC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLSC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATCONT.storedata=1;

OCMATCONT.HE.equilibriumcoord=1:length(limSetdepvar);
switch pathtpe
    case 's'
        OCMATLSC.subspacedim=OCMATLSC.numstable;
    case 'u'
        OCMATLSC.subspacedim=OCMATLSC.numunstable;
    case {'sc','cs'}
        OCMATLSC.subspacedim=OCMATLSC.numstable+OCMATLSC.numcenter;
end
OCMATLSC.orthspacedim=limSetDim-OCMATLSC.subspacedim;
Y=zeros(OCMATLSC.orthspacedim,OCMATLSC.subspacedim);
OCMATLSC.Y=Y;
% orthogonal basis for stable eigenspace
OCMATLSC.Q0=computeBase(JacobianMatrix,stableflag,OCMATLSC.subspacedim);
OCMATCONT.HE.Ycoord=reshape(length(limSetdepvar)+(1:OCMATLSC.orthspacedim*OCMATLSC.subspacedim),OCMATLSC.orthspacedim,OCMATLSC.subspacedim);
OCMATLSC.Id=eye(OCMATLSC.orthspacedim);
OCMATLSC.numY=numel(Y);

sol.parameters=[sol.parameters limSetdepvar.' Y(:).'];

if freeendtime
    OCMATLSC.freeendtimecoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
    if freeendtime>0
        depvar=dependentvar(ocAsym);
        OCMATLSC.distance=norm(OCMATLSC.saddlepoint-depvar(:,end));
    end
else
    OCMATLSC.truncationtimecoord=[];
end
OCMATLSC.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
sol.parameters=[sol.parameters OCMATLSC.parametervalue(parindex)];
OCMATLSC.parametercoord=length(sol.parameters);

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=size(OCMATLSC.asymptoticmatrix,2);

OCMATLSC.objectivevaluecalc=objectivevaluecalc;
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocAsym,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocAsym,1)))];
end
if objectivevaluecalc
    OCMATLSC.objectivevaluecoord=size(sol.y,1);
end

%compute borders
OCMATCONT.OPTIONS.xyvectorized=1;
OCMATCONT.bvpmethod=ocAsym.solver;


ocmatmsg('\nTo initialize the bordered matrix set\n\topt=setocoptions(''OCCONTARG'',''WorkSpace'',1);\n')
w=[];
v=[];

OCMATLSC.LAE_phi=w(:);
OCMATLSC.LAE_psi=v(:);
OCMATLSC.LAE_switch=1;

OCMATLSC.freeendtime=freeendtime;
OCMATLSC.autonomous=isautonomous(ocObj);
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



function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
sol.solver=solvername;
