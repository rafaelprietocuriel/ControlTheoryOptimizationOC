function sol=initocmat_LC_H_PII(ocObj,ocEP,parindex,amp,opt,varargin)
% INITOCMAT_LC_H_P initialization for the continuation of a limit cycle
% with respect to a parameter, starting from a Hopf bifurcation
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX)
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATLC
global OCMATCONT OCMATLC
sol=[];
targetparametervalue=[];
linearizationcalc=[];
fixcoordinate=[];
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
monodromyidx=find(strcmpi(varargin,'monodromy'));
fixcoordinateidx=find(strcmpi(varargin,'fixcoordinate'));

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocEP)
    ocmatmsg('oc equilibrium is empty.')
    return
end
if ~isequilibrium(ocEP)
    ocmatmsg('Second argument is mot an equilibrium.')
    return
end
if nargin==4
    opt=[];
end
if isempty(opt)
    opt=defaultocoptions;
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(fixcoordinateidx)
    fixcoordinate=varargin{fixcoordinateidx+1};
end
if isempty(linearizationcalc)
    linearizationcalc=0;
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(targetparametervalue) && length(parindex)>1
    targetparametervalue=zeros(length(parindex),1);
end
if ~isempty(monodromyidx)
    monodromy=varargin{monodromyidx+1};
else
    monodromy=0;
end

targetparametervalue=targetparametervalue(:).';

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modeltype=modeltype(ocObj);
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4LimitCycleContinuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

switch OCMATCONT.modeltype
    case 'standardmodel'
        OCMATLC.canonicalsystem=funch{1};
        OCMATLC.canonicalsystemjacobian=funch{2}{1};
        OCMATLC.canonicalsystemparameterjacobian=funch{2}{2};
        OCMATLC.canonicalsystemhessian=funch{3}{1};
        OCMATLC.canonicalsystemparameterhessian=funch{3}{2};
        % function for the boundary conditions
        OCMATLC.bcperiodic=funch{5}{1};

        % function for Jacobian
        OCMATLC.bcjacobianperiodic=funch{6}{1};

        % function describing the hybrid structure of the problem
        OCMATLC.hybridinfo=funch{7}{1};
        OCMATLC.domain=funch{7}{2};
        OCMATLC.guard=funch{7}{3};
        OCMATLC.reset=funch{7}{4};
        OCMATLC.switchtime=funch{7}{5};
        OCMATLC.jacobianguard=funch{7}{7};
        OCMATLC.jacobianreset=funch{7}{8};
        OCMATLC.domaindiscretization=funch{7}{9};
        OCMATLC.timesettransformation=funch{7}{10};
        if ~isautonomous(ocObj)
            OCMATLC.canonicalsystemderivativetime=funch{2}{3};
        end
        hybridinfo=OCMATLC.hybridinfo();
        % retrieve information about the arcs from the specifc model file
        % domaindata and discretizationdata are ordered in the sequence the
        % arcarguments appear in hybridinfo.arcarg
        for ii=1:numel(hybridinfo.arcarg)
            domaindata(ii)=OCMATLC.domain(hybridinfo.arcarg(ii));
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
    case 'odemodel'

        OCMATLC.dynamics=funch{1};
        OCMATLC.dynamicsjacobian=funch{2}{1};
        OCMATLC.dynamicsparameterjacobian=funch{2}{2};
        OCMATLC.dynamicshessian=funch{3}{1};
        OCMATLC.dynamicsparameterhessian=funch{3}{2};
        % function for the boundary conditions
        OCMATLC.bcperiodic=funch{5}{1};

        % function for Jacobian
        OCMATLC.bcjacobianperiodic=funch{6}{1};
end
% general function
OCMATLC.plotcontinuation=funch{11};
OCMATLC.testadmissibility=funch{12};
OCMATLC.datapath=funch{20};
OCMATLC.saveintermediatefiles=funch{21};
OCMATLC.fixstep=0;


J=jacobian(ocEP);
% calculate eigenvalues and eigenvectors
[V,D] = eig(J);
% find pair of complex eigenvalues
d = diag(D);
idx=find(imag(d)==0);

smallest_sum = Inf;
numode=size(J,1);
searchidx=setdiff(1:numode-1,idx);
for ii=searchidx(1:end-1)
  [val,idx] = min(abs(d(ii+1:numode)+d(ii)));
  if val < smallest_sum
    idx1 = ii;
    idx2 = ii+idx;
    smallest_sum = val;
  end
end
% real part? Oh dear, a neutral saddle!
if imag(d(idx1)) == 0 && imag(d(idx2)) == 0
  ocLC=dynprimitive([]);
  ocmatmsg('Neutral saddle\n');
  return;
end
% get imaginary part and corresponding eigenvector
omega = abs(imag(d(idx2)));
Q = V(:,idx1);

d = real(Q)'*real(Q);
s = imag(Q)'*imag(Q);
r = real(Q)'*imag(Q);
Q = Q*exp(i*atan2(2*r,s-d)/2);
Q = Q/norm(real(Q));

% initial amplitude h
% calculate initial cycle and its tangent vector
lc.octrajectory.x=linspace(0,1,opt.GENERAL.TrivialArcMeshNum);
oldlc.x=linspace(0,1,opt.GENERAL.TrivialArcMeshNum);
y=kron(exp(2*pi*i*linspace(0,1,opt.GENERAL.TrivialArcMeshNum)),Q);
y=reshape(y,numode,[]);
oldlc.yp=-imag(y);
oldlc.solver='';
OCMATLC.oldlcInit=oldlc;
OCMATLC.oldlc=oldlc;
tangent = real(y);
lc.octrajectory.y=repmat(ocEP.y,1,opt.GENERAL.TrivialArcMeshNum)+amp*tangent;
lc.octrajectory.arcarg=arcargument(ocEP);
lc.octrajectory.arcposition=[1;opt.GENERAL.TrivialArcMeshNum];
lc.octrajectory.arcinterval=[0 2*pi/omega];
lc.octrajectory.solverinfo.tmesh=[lc.octrajectory.x];
if length(parindex)==1 && isempty(targetparametervalue)
    lc.octrajectory.solverinfo.coeff=[lc.octrajectory.y(:);2*pi/omega;parametervalue(ocObj,parindex).'];
    lc.octrajectory.solverinfo.parameters=[2*pi/omega;parametervalue(ocObj,parindex).'];
    lc.octrajectory.solverinfo.tangent=[tangent(:);0;0];
else
    lc.octrajectory.solverinfo.coeff=[lc.octrajectory.y(:);2*pi/omega;0];
    lc.octrajectory.solverinfo.parameters=[2*pi/omega;0];
    lc.octrajectory.solverinfo.tangent=[tangent(:);1];
end
lc.octrajectory.solverinfo.tangent=lc.octrajectory.solverinfo.tangent/norm(lc.octrajectory.solverinfo.tangent);
lc.period=2*pi/omega;
ocLC=dynprimitive(lc);

% mode and path specific variables
depvar=dependentvar(ocLC);
OCMATLC.parametervalue=parametervalue(ocObj);
OCMATLC.autonomous=isautonomous(ocObj);
OCMATLC.initialtime=0;
OCMATLC.switchtimecoord=[];
OCMATLC.varyparameterindex=parindex;
OCMATLC.periodcoord=1;
OCMATLC.parametercoord=2;
OCMATLC.targetparametervalue=targetparametervalue;
if ~isempty(OCMATLC.targetparametervalue)
    OCMATLC.initialparametervalue=OCMATLC.parametervalue(parindex);
    OCMATLC.initialparametervalue=OCMATLC.initialparametervalue(:).';
    OCMATLC.continuationvector=OCMATLC.targetparametervalue-OCMATLC.initialparametervalue;
end
OCMATLC.linearizationcalc=linearizationcalc;
OCMATLC.fixcoordinate=fixcoordinate;
if isempty(fixcoordinate)
    OCMATLC.velocityvector=oldlc.yp(:,1);
    OCMATLC.velocitycoord=1:length(oldlc.yp(:,1));
    OCMATLC.velocityvector=OCMATLC.velocityvector/norm(OCMATLC.velocityvector);
else
    OCMATLC.velocityvector=[];
    OCMATLC.fixcoordinate=fixcoordinate;
    OCMATLC.fixvalue=depvar(fixcoordinate,1);
end
OCMATLC.initialpoint=depvar(:,1);
OCMATLC.monodromy=monodromy;
OCMATLC.objectivevaluecalc=0;

OCMATCONT.codimension=1;
OCMATCONT.numintegralconstraint=0;

pathname=OCMATLC.datapath();
[resultfile,globalvarfile]=OCMATLC.saveintermediatefiles();
OCMATLC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
sol=generatesolstruct(ocLC);

if ~strcmp(OCMATCONT.modeltype,'standardmodel')
    return
end
% test if pure state constraints are defined
OCMATLC.stateconstraint=stateconstraint(ocObj);


function sol=generatesolstruct(ocLC,solvername,varargin)

sol.x=independentvar(ocLC);
sol.y=dependentvar(ocLC);
sol.arcarg=arcargument(ocLC);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocLC);
sol.arcposition=arcposition(ocLC);
sol.parameters=parameters(ocLC);
sol.solver='';

sol.solverinfo.tangent=tangent(ocLC);
sol.solverinfo.coeff=ocLC.solverinfo.coeff;
sol.solverinfo.tmesh=ocLC.solverinfo.tmesh;
