function sol=initocmat_DAE_DAE(docObj,docAsym,contcoordinate,targetvalue,opt)
%
% INITOCMAT_DAE_DAE initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_DAE_DAE(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_DAE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_DAE_EP
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtpe=pathtype(docAsym);
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
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(docObj);
OCMATCONT.modelfunc=modelspecificfunc(docObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
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
J=explicitjacobian(limSet);
OCMATAE.parametervalue=parametervalue(docObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.linearization=linearization(limSet);
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.saddlepoint=OCMATAE.saddlepoint(:,1);
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,'d',ZeroDeviationTolerance);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=sol.y(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATCONT.numarc=length(sol.arcarg);


function sol=generateodiffestruct(docAsym)

sol.x=[docAsym.x0 docAsym.x];
sol.y=[docAsym.y0 docAsym.y];
sol.x0=docAsym.x0;
% add continuation parameter value
sol.parameters=[docAsym.solverinfo.parameters(1:end-1) 0];
sol.arcarg=arcargument(docAsym);
sol.arcposition=docAsym.arcposition;
sol.arcposition(2,end)=sol.arcposition(2,end)+1;