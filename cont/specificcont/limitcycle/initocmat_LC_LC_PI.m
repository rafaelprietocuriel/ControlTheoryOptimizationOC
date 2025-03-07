function sol=initocmat_LC_LC_PI(ocObj,ocLC,parindex,varargin)
% INITOCMAT_LC_H_P initialization for the continuation of a limit cycle
% with respect to a parameter, starting from a Hopf bifurcation
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX)
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATLC
global OCMATCONT OCMATLC
sol=[];
monodromy=[]; % provide possible jumps of costate for pure state constraints
targetparametervalue=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocLC)
    ocmatmsg('oc limitcycle is empty.')
    return
end
if ~isperiodic(ocLC)
    ocmatmsg('Second argument is not a limitcycle.')
    return
end
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
monodromyidx=find(strcmpi(varargin,'monodromy'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(monodromyidx)
    monodromy=varargin{monodromyidx+1};
end

if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4LimitCycleContinuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

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

sol=generatesolstruct(ocLC);
numswitchtimes=length(sol.arcinterval(2:end-1));

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    if monodromy
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+OCMATCONT.DOMAINDDATA(ii).numeq^2;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
sol.parameters=[sol.parameters parametervalue(ocObj,parindex)];

% mode and path specific variables
depvar=dependentvar(ocLC);
if isocmodel(ocObj)
    dxdt=canonicalsystem(ocObj,ocLC,[],1);
elseif isodemodel(ocObj)
    dxdt=dynamics(ocObj,ocLC,[],1);
end
OCMATLC.parametervalue=parametervalue(ocObj);
OCMATLC.autonomous=isautonomous(ocObj);
OCMATLC.initialtime=sol.x0;
OCMATLC.switchtimecoord=1:numswitchtimes;
OCMATLC.periodcoord=1+numswitchtimes;
OCMATLC.parametercoord=2+numswitchtimes;
OCMATLC.velocityvector=dxdt(:,1);
OCMATLC.velocitycoord=1:length(dxdt(:,1));
OCMATLC.velocityvector=OCMATLC.velocityvector/norm(OCMATLC.velocityvector);
OCMATLC.initialpoint=depvar(:,1);
OCMATLC.varyparameterindex=parindex;
OCMATLC.targetparametervalue=targetparametervalue;
OCMATLC.calcmonodromy=monodromy;

n=length(sol.y(:,1));
OCMATLC.depvarcoord=1:n;
OCMATLC.numdepvar=n;
if monodromy
    id=eye(length(depvar(:,1)));
    sol.y=[sol.y;repmat(id(:),1,length(sol.x))];
    OCMATLC.id=id(:);
    OCMATLC.monodromycoord=n+1:length(sol.y(:,1));
end
OCMATCONT.codimension=1;
OCMATCONT.numintegralconstraint=0;

pathname=OCMATLC.datapath();
[resultfile,globalvarfile]=OCMATLC.saveintermediatefiles();
OCMATLC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

function sol=generatesolstruct(ocLC,varargin)

sol.x=independentvar(ocLC);
sol.y=dependentvar(ocLC);
sol.arcarg=arcargument(ocLC);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocLC);
sol.arcposition=arcposition(ocLC);
sol.parameters=parameters(ocLC);
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end);
else
    numcontpar=length(continuationparameter(ocLC));
    sol.parameters(end-numcontpar+1:end)=[];

end
sol.solver='';

try
sol.solverinfo.tangent=tangent(ocLC);
sol.solverinfo.coeff=ocLC.solverinfo.coeff;
sol.solverinfo.tmesh=ocLC.solverinfo.tmesh;
end