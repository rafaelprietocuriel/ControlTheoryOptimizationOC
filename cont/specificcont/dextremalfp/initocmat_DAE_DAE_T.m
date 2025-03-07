function sol=initocmat_DAE_DAE_T(ocObj,ocAsym,timeindex,initialcoordinate,opt,varargin)
% INITOCMAT_DAE_EP initialization for asymptotic extremal calculation
%
% ocObj ...   corresponding optimal control model
% ocAsym ...  manifold path (ocasymptotic)
% contidx ... coordinates of the continuation variable
% targetvalue ... determines direction of the continuation
%
% SOL=INITOCMAT_DAE_DAE(ocObj,ocAsym,contidx,targetvalue,discretizationdata)
% discretizationdata specifies discretization data (necessary fields depend on the
% solver)
%       discretizationdata.solver ... the corresponding bvp solver
%       discretizationdata.nummesh   ... number of mesh points
%       discretizationdata.asymptoticapproximation ... real number means truncation of
%               the integration time, 'inf' transformation to unit interval
%       discretizationdata.normalization ... 1 integration interval normalized to one, 0
%               no normalization of the integration interval
%       discretizationdata.numcols ... number of collocation
%                      points
%       discretizationdata.colmethod ... 'gaussian', 'lobatto'
clear global OCMATCONT OCMATDAE
global OCMATCONT OCMATDAE
sol=[];
pathtpe=pathtype(ocAsym);
targetvalue=[];
targettype='';
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
targettypeidx=find(strcmpi(varargin,'targettype'));
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if isempty(targetvalue)
    targetvalue=1e-1;
end
if isempty(targettype)
    targettype='d';
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATDAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATDAE.canonicalsystem=funch{1};
OCMATDAE.canonicalsystemjacobian=funch{2}{1};
OCMATDAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATDAE.canonicalsystemhessian=funch{3}{1};
OCMATDAE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATDAE.bcinitial=funch{5}{1};
OCMATDAE.bcasymptotic=funch{5}{2};
OCMATDAE.bctransversality=funch{5}{3};
OCMATDAE.equilibrium=funch{5}{4};

% function for Jacobian
OCMATDAE.bcjacobianinitial=funch{6}{1};
OCMATDAE.bcjacobianasymptotic=funch{6}{2};
OCMATDAE.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATDAE.hybridinfo=funch{7}{1};
OCMATDAE.domain=funch{7}{2};
OCMATDAE.guard=funch{7}{3};
OCMATDAE.reset=funch{7}{4};
OCMATDAE.switchtime=funch{7}{5};
OCMATDAE.jacobianguard=funch{7}{7};
OCMATDAE.jacobianreset=funch{7}{8};
OCMATDAE.domaindiscretization=funch{7}{9};
OCMATDAE.timesettransformation=funch{7}{10};

% general function
OCMATDAE.plotcontinuation=funch{11};
OCMATDAE.testadmissibility=funch{12};
OCMATDAE.datapath=funch{20};
OCMATDAE.saveintermediatefiles=funch{21};

hybridinfo=OCMATDAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATDAE.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocAsym,solver(ocAsym));
arctimes=arcinterval(ocAsym);
sol.parameters=[sol.parameters arctimes(timeindex)];
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
limSet=limitset(ocAsym);
depvar=dependentvar(ocAsym);
OCMATDAE.parametervalue=parametervalue(ocObj);
OCMATDAE.initialtime=arctimes(1);
OCMATDAE.switchtimecoord=1:numel(sol.parameters)-1;
OCMATDAE.inftimetransformation=inftimetransformation(ocAsym);
OCMATDAE.truncationtime=arctimes(end);
OCMATDAE.linearization=linearization(limSet);
OCMATDAE.saddlepoint=dependentvar(limSet);
OCMATDAE.asymptoticmatrix=asymptoticbc(ocAsym.limitset.linearization,pathtpe,'c',ZeroDeviationTolerance);
OCMATDAE.pathtype=pathtpe;
pathname=OCMATDAE.datapath();
[resultfile,globalvarfile]=OCMATDAE.saveintermediatefiles();
OCMATDAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATDAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATDAE.initialstate=depvar(initialcoordinate,1);
OCMATDAE.initialcoordinate=initialcoordinate;
OCMATDAE.varyparameterindex=timeindex;
OCMATDAE.targetvalue=targetvalue;
OCMATDAE.targettype=targettype;

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATDAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;