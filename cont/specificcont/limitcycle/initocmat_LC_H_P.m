function sol=initocmat_LC_H_P(ocObj,ocEP,parindex,amp,opt,varargin)
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
if nargin>=6
    targetparametervalue=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
linearizationidx=find(strcmpi(varargin,'linearization'));
if ~isempty(linearizationidx)
    linearizationcalc=varargin{linearizationidx+1};
end
if isempty(linearizationcalc)
    linearizationcalc=0;
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end

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

hybridinfo=OCMATLC.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATLC.domain(hybridinfo.arcarg(ii));
end

J=jacobian(ocEP);
% calculate eigenvalues and eigenvectors
[V,D] = eig(J);
% find pair of complex eigenvalues
d = diag(D);
smallest_sum = Inf;
numode=size(J,1);
for ii=1:numode-1
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
%oldlc.yp=-imag(y);
oldlc.solver='';
tangent = real(y);
lc.octrajectory.y=repmat(ocEP.y,1,opt.GENERAL.TrivialArcMeshNum)+amp*tangent;
lc.octrajectory.arcarg=arcargument(ocEP);
lc.octrajectory.arcposition=[1;opt.GENERAL.TrivialArcMeshNum];
lc.octrajectory.arcinterval=[0 2*pi/omega];
lc.octrajectory.solverinfo.tmesh=[lc.octrajectory.x];
lc.octrajectory.solverinfo.coeff=[lc.octrajectory.y(:);2*pi/omega;parametervalue(ocObj,parindex)];
lc.octrajectory.solverinfo.tangent=[tangent(:);0;0];
lc.octrajectory.solverinfo.tangent=lc.octrajectory.solverinfo.tangent/norm(lc.octrajectory.solverinfo.tangent);
lc.octrajectory.solverinfo.parameters=[2*pi/omega;parametervalue(ocObj,parindex)];
lc.period=2*pi/omega;
ocLC=dynprimitive(lc);
oldlc.yp=2*pi/omega*canonicalsystem(ocObj,ocLC);
OCMATLC.oldlcInit=oldlc;
OCMATLC.oldlc=oldlc;

% test if pure state constraints are defined
OCMATLC.stateconstraint=stateconstraint(ocObj);
sol=generatesolstruct(ocLC);
sccounter=0;
if OCMATLC.stateconstraint
    arcarg4sc=arcargwithactivestateconstraint(ocObj);
    for ii=2:length(sol.arcarg)
        if ismember(sol.arcarg(ii),arcarg4sc)
            sccounter=sccounter+1;
            OCMATLC.jumpcostateindex(sccounter)=ii;
        end
    end
    if sccounter
        OCMATLC.jumpcostatecoord=1:sccounter;
    else
        OCMATLC.jumpcostatecoord=[];
        OCMATLC.jumpcostateindex=[];
    end
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

% mode and path specific variables
OCMATLC.parametervalue=parametervalue(ocObj);
OCMATLC.autonomous=isautonomous(ocObj);
OCMATLC.initialtime=sol.x0;
OCMATLC.switchtimecoord=[];
OCMATLC.varyparameterindex=parindex;
OCMATLC.periodcoord=1;
OCMATLC.parametercoord=2;
OCMATLC.targetparametervalue=targetparametervalue;
OCMATLC.linearizationcalc=linearizationcalc;
if linearizationcalc
    OCMATLC.linearizationinit=eye(2*statenum(ocObj));
    OCMATLC.linearizationinit=OCMATLC.linearizationinit(:);
    OCMATLC.linearizationindex=2*statenum(ocObj)+[1:4*statenum(ocObj)^2];
    OCMATLC.linearizationindex=reshape(OCMATLC.linearizationindex,2*statenum(ocObj),2*statenum(ocObj));
end
OCMATCONT.codimension=1;
OCMATCONT.numintegralconstraint=1;

pathname=OCMATLC.datapath();
[resultfile,globalvarfile]=OCMATLC.saveintermediatefiles();
OCMATLC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

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
