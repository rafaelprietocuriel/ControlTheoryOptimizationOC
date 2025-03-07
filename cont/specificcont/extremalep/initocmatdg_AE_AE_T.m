function sol=initocmatdg_AE_AE_T(dgObj,ocAsym,timeindex,initialcoordinate,opt,varargin)
% INITOCMAT_AE_EP initialization for asymptotic extremal calculation
%
% dgObj ...   corresponding optimal control model
% ocAsym ...  manifold path (ocasymptotic)
% contidx ... coordinates of the continuation variable
% targetvalue ... determines direction of the continuation
%
% SOL=INITOCMAT_AE_AE(dgObj,ocAsym,contidx,targetvalue,discretizationdata)
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
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtpe=pathtype(ocAsym);
targetvalue=[];
targettype='';
objectivevaluecalc=[];
if isempty(dgObj)
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
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));

if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end

if isempty(targetvalue)
    targetvalue=1e-1;
end
if isempty(targettype)
    targettype='T';
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(dgObj);
OCMATCONT.modelfunc=modelspecificfunc(dgObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
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

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

sol=generatesolstruct(ocAsym,solver(ocAsym));
arctimes=arcinterval(ocAsym);
if isempty(timeindex)
    timeindex=length(arctimes);
end

numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
end

sol.parameters=[sol.parameters arctimes(timeindex)];
OCMATAE.autonomous=isautonomous(dgObj);


% mode and path specific variables
limSet=limitset(ocAsym);
depvar=dependentvar(ocAsym);
J=linearization(limSet);
OCMATAE.linearization=J;

OCMATAE.parametervalue=parametervalue(dgObj);
OCMATAE.initialtime=arctimes(1);
OCMATAE.switchtimecoord=1:numswitchtimes;
OCMATAE.truncationtime=arctimes(end);
OCMATAE.saddlepoint=dependentvar(limSet);
%OCMATAE.asymptoticmatrix=asymptoticbc(OCMATAE.linearization,pathtpe,'c',ZeroDeviationTolerance);
if isequilibrium(limSet)
    limtype='c';
else
    limtype='d';
end
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,limtype,ZeroDeviationTolerance);
OCMATAE.objectivevaluecoord=objectivevaluecoordinate(ocAsym);
if objectivevaluecalc && isempty(OCMATAE.objectivevaluecoord)
    o=objectivefunction(dgObj,ocAsym,1);
    OCMATAE.objectivevaluecoord=size(sol.y,1)+(1:OCMATAE.playernum);
    sol.y(end+(1:OCMATAE.playernum),:)=[zeros(OCMATAE.playernum,1) cumsum((o(:,1:end-1)+o(:,2:end))/2.*repmat(diff(time(dgObj,sol,1)),OCMATAE.playernum,1))];
end

OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=objectivevaluecalc;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATAE.initialstate=depvar(initialcoordinate,1);
OCMATAE.initialcoordinate=initialcoordinate;
OCMATAE.varyparameterindex=timeindex;
OCMATAE.targetvalue=targetvalue;
OCMATAE.targettype=targettype;
OCMATAE.movinghorizon=1;
OCMATAE.movinghorizoncoord=length(sol.parameters);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATCONT.numeq=size(sol.y,1);
OCMATCONT.canonicalsystemcoordinate=1:canonicalsystemdimension(dgObj);

OCMATCONT.equilibriumcoord=1:length(dependentvar(limitset(ocAsym)));


function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
sol.solver=solvername;

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