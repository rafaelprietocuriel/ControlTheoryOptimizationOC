function sol=initocmat_DAE_DAE_P(docObj,docAsym,parindex,initialcoordinate,opt,varargin)
% initocmat_DAE_DAE_P initialization for asymptotic extremal calculation
% varying a parameter
%
% SOL=initocmat_DAE_DAE_P(docObj,ocEP,contidx,targetvalue)


clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtpe=pathtype(docAsym);
targetparametervalue=[];
if isempty(docObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(docAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
if nargin>=6
    targetparametervalue=varargin{1};
end
stableflag=strcmp(pathtpe,'s');

if ischar(parindex)
    parindex=parameterindex(docObj,parindex);
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(docObj);
OCMATCONT.modelfunc=modelspecificfunc(docObj,'4SaddlePathContinuation');

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(docObj);
OCMATCONT.modelfunc=modelspecificfunc(docObj,'4SaddlePathContinuation');
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystemmap=funch{1};
OCMATAE.canonicalsystemmapjacobian=funch{2}{1};
OCMATAE.canonicalsystemmapparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemmaphessian=funch{3}{1};
OCMATAE.canonicalsystemmapparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.equilibrium=funch{5}{4};

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianasymptotic=funch{6}{2};
OCMATAE.bcjacobiantransversality=funch{6}{3};

% general function
OCMATAE.findarcposition=funch{10};
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

sol=generateodiffestruct(docAsym);

% mode and path specific variables
limSet=limitset(docAsym);
JacobianMatrix=explicitjacobian(limSet);

OCMATAE.parametervalue=parametervalue(docObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.linearization=linearization(limSet);
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.pathtype=pathtpe;
OCMATAE.stableflag=stableflag;
[OCMATAE.asymptoticmatrix OCMATAE.numstable OCMATAE.numunstable OCMATAE.numcenter infoStruct]=asymptoticbc(JacobianMatrix,pathtpe,'d',ZeroDeviationTolerance);
OCMATAE.initialstate=sol.y(initialcoordinate,1);
OCMATAE.initialcoordinate=initialcoordinate;
OCMATAE.varyparameterindex=parindex;
OCMATAE.targetparametervalue=targetparametervalue;

OCMATCONT.codimension=1;

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

Y=zeros(OCMATAE.numunstable,OCMATAE.numstable);
OCMATAE.Y=Y;
OCMATAE.Q0=computeBase(JacobianMatrix,stableflag,OCMATAE.numstable);

sol.parameters=[sol.parameters dependentvar(limSet).' Y(:).' OCMATAE.parametervalue(parindex)];
sol.solverinfo.coeff=[sol.y(:);dependentvar(limSet);Y(:);OCMATAE.parametervalue(parindex)];
sol.solverinfo.tangent=[];%tangent(docAsym);

OCMATCONT.HE.equilibriumcoord=1:numel(dependentvar(limSet));
if stableflag
    dimSubSpace=OCMATAE.numstable;
    OCMATCONT.HE.Ycoord=reshape(numel(dependentvar(limSet))+(1:OCMATAE.numstable*OCMATAE.numunstable),OCMATAE.numstable,OCMATAE.numunstable);
    OCMATAE.Id=eye(OCMATAE.numunstable);
else
    dimSubSpace=OCMATAE.numunstable;
    OCMATCONT.HE.Ycoord=reshape(numel(dependentvar(limSet))+(1:OCMATAE.numstable*OCMATAE.numunstable),OCMATAE.numunstable,OCMATAE.numstable);
    OCMATAE.Id=eye(OCMATAE.numstable);
end
OCMATAE.dimSubSpace=dimSubSpace;
OCMATAE.numY=numel(Y);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);
OCMATCONT.HE.Ycoord=OCMATCONT.HE.Ycoord+numel(sol.arcarg)-1;
OCMATCONT.HE.equilibriumcoord=OCMATCONT.HE.equilibriumcoord+numel(sol.arcarg)-1;
OCMATCONT.numarc=length(sol.arcarg);

function  Q=computeBase(J,stableflag,NSub)

if ~stableflag
    J=-J;
end
[VU, DU] = realeig(J);
% Select first NSub eigenvectors: unstable eigenspace
% Compute orthonormal basis for the eigenspace
VU = VU(:,1:NSub);
[Q,RU] = qr(VU);



function sol=generateodiffestruct(docAsym)

sol.x=[docAsym.x0 docAsym.x];
sol.y=[docAsym.y0 docAsym.y];
sol.x0=docAsym.x0;
% add continuation parameter value
sol.parameters=[docAsym.solverinfo.parameters(1:end-1)];
sol.arcarg=arcargument(docAsym);
sol.arcposition=docAsym.arcposition;
sol.arcposition(2,end)=sol.arcposition(2,end)+1;