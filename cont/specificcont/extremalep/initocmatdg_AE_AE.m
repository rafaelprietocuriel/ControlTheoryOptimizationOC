 function sol=initocmatdg_AE_AE(dgObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
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
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
movinghorizon=[];
pathtpe=pathtype(ocAsym);
objectivevaluecalc=[];
targetstate=[];
targetstatecoordinate=[];
if isempty(pathtpe)
    pathtpe='s';
end
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
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
pathtypeidx=find(strcmpi(varargin,'pathtype'));
targetstateidx=find(strcmpi(varargin,'targetstate'));
targetstatecoordinateidx=find(strcmpi(varargin,'targetstatecoordinate'));
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(targetstateidx)
    targetstate=varargin{targetstateidx+1};
end
if ~isempty(targetstatecoordinateidx)
    targetstatecoordinate=varargin{targetstatecoordinateidx+1};
end
if isempty(movinghorizon)
    movinghorizon=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if ~isempty(pathtypeidx)
    pathtpe=varargin{pathtypeidx+1};
end
targetvalue=targetvalue(:);
OCMATAE.targetstate=targetstate;
OCMATAE.targetstatecoordinate=targetstatecoordinate;

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
OCMATAE.freeendtime=[];

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(dgObj);
OCMATCONT.modelfunc=modelspecificfunc(dgObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
if ~isautonomous(dgObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
end
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
% if stopcriterion
%     OCMATAE.stopcriterionfunc=funch{5}{6};
% end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

OCMATAE.playernum=playernum(dgObj);

sol=generatesolstruct(ocAsym,solver(ocAsym));

numswitchtimes=length(sol.arcinterval(2:end-1));
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
elseif length(sol.parameters)==numswitchtimes
    sol.parameters=[sol.parameters];
elseif length(sol.parameters)<numswitchtimes
    sol.parameters=sol.arcinterval(2:end-1);
end

% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
OCMATAE.linearization=J;

depvar=dependentvar(ocAsym);
OCMATAE.parametervalue=parametervalue(dgObj);
OCMATAE.autonomous=isautonomous(dgObj);

OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numswitchtimes;
switch pathtpe
    case {'s','sc','cs','sts','ws'}
        OCMATAE.truncationtime=sol.arcinterval(end);
    case {'u','uc','cu','wu','stu'}
        OCMATAE.truncationtime=-abs(sol.arcinterval(end));
        sol.arcinterval=-abs(sol.arcinterval(end));
end
OCMATAE.saddlepoint=dependentvar(limSet);
asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.asymptoticmatrix=asymptoticmatrix;

OCMATCONT.equilibriumcoord=1:length(dependentvar(limitset(ocAsym)));

OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
OCMATAE.objectivevaluecalc=objectivevaluecalc;
OCMATAE.movinghorizon=movinghorizon;
continuationparameter=0;
if movinghorizon
    yend=depvar(OCMATCONT.equilibriumcoord,end);
    hatx=OCMATAE.saddlepoint;
    OCMATAE.distance=norm(yend-hatx);
    OCMATAE.movinghorizoncoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
end
sol.parameters=[sol.parameters continuationparameter];
OCMATAE.objectivevaluecoord=objectivevaluecoordinate(ocAsym);
if objectivevaluecalc && isempty(OCMATAE.objectivevaluecoord)
    o=objectivefunction(dgObj,ocAsym,1);
    OCMATAE.objectivevaluecoord=size(sol.y,1)+(1:OCMATAE.playernum);
    sol.y(end+(1:OCMATAE.playernum),:)=[zeros(OCMATAE.playernum,1) cumsum((o(:,1:end-1)+o(:,2:end))/2.*repmat(diff(time(dgObj,sol,1)),OCMATAE.playernum,1))];
end

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;
OCMATCONT.numeq=size(sol.y,1);
OCMATCONT.canonicalsystemcoordinate=1:canonicalsystemdimension(dgObj);

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
if isfield(ocTrj.solverinfo,'conttype')
    switch ocTrj.solverinfo.conttype
        case 'extremal2ep'
            if ~isempty(sol.parameters)
                numcontpar=length(continuationparameter(ocTrj));
                sol.parameters(end-numcontpar+1:end)=[];
            end
        otherwise
            sol.parameters=[];
    end
end
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