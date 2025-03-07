function sol=initocmat_HOM_LCI(ocObj,ocHom,parindex,opt,varargin)
% INITOCMAT_LC_H_P initialization for the continuation of a limit cycle
% with respect to a parameter, starting from a Hopf bifurcation
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX)
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATHET
global OCMATCONT OCMATHET
sol=[];

if nargin==3
    opt=defaultocoptions;
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocHom)
    ocmatmsg('oc limitcycle is empty.')
    return
end
if nargin==4
    opt=[];
end
if isempty(opt)
    opt=defaultocoptions;
end
ocHom=ocmultipath(ocHom);
hetorder=multiplicity(ocHom);
for ii=1:2
    pathtpe{ii}=pathtype(ocHom(ii));
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATHET.canonicalsystem=funch{1};
OCMATHET.canonicalsystemjacobian=funch{2}{1};
OCMATHET.canonicalsystemparameterjacobian=funch{2}{2};
OCMATHET.canonicalsystemhessian=funch{3}{1};
OCMATHET.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions

% function for Jacobian

% function describing the hybrid structure of the problem
OCMATHET.hybridinfo=funch{7}{1};
OCMATHET.domain=funch{7}{2};
OCMATHET.guard=funch{7}{3};
OCMATHET.reset=funch{7}{4};
OCMATHET.switchtime=funch{7}{5};
OCMATHET.jacobianguard=funch{7}{7};
OCMATHET.jacobianreset=funch{7}{8};
OCMATHET.domaindiscretization=funch{7}{9};
OCMATHET.timesettransformation=funch{7}{10};
if ~isautonomous(ocObj)
    OCMATHET.canonicalsystemderivativetime=funch{2}{3};
end

% function for the boundary conditions
OCMATHET.bcasymptotic=funch{5}{2};
% general function
OCMATHET.plotcontinuation=funch{11};
OCMATHET.testadmissibility=funch{12};
OCMATHET.datapath=funch{20};
OCMATHET.saveintermediatefiles=funch{21};

hybridinfo=OCMATHET.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATHET.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocHom);
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
OCMATHET.hetorder=hetorder;
OCMATCONT.HE.edge=[];
arcoffset=0;
switchtimes=[];
for ii=1:2
    arcn=arcnum(ocHom(ii));
    %solverInfoStruct=solverinfo(ocHom(ii));
    OCMATHET.switchtimecoord{ii}=length(switchtimes)+(1:arcn-1);
    switchtimes=[switchtimes ocHom(ii).arcinterval(2:end-1)];
    arcarg=arcargument(ocHom(ii));

    OCMATHET.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATHET.numarc(ii)=arcn;
    OCMATHET.initialstateindex(ii)=numel(ocHom(ii).x);
    OCMATHET.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcn;
end
OCMATHET.initialstateindex=cumsum(OCMATHET.initialstateindex);
OCMATHET.initialstateindex=[0 OCMATHET.initialstateindex(end-1)]+1;
OCMATHET.solutionindex=zeros(1,sum(OCMATHET.numarc));% relate arcindex to indifference solution
counter=0;
for ii=1:2
    counter_start=counter+1;
    counter=counter+OCMATHET.numarc(ii);
    OCMATHET.solutionindex(counter_start:counter)=ii;
end

OCMATHET.cumsumnumarc=cumsum(OCMATHET.numarc);
OCMATHET.initcoord=[1 OCMATHET.cumsumnumarc(1:end-1)+1];

OCMATHET.parametervalue=parametervalue(ocObj);
OCMATHET.initialtime=sol.x0;

righttimeindex=cumsum(OCMATHET.numarc+1);
OCMATHET.truncationtime=sol.arcinterval(righttimeindex);

% phase condition
dxdt=canonicalsystem(ocObj,ocHom(1),[],1);
depvar=dependentvar(ocHom(1));
depvar2=dependentvar(ocHom(2));
OCMATHET.velocityvector=dxdt(:,1);
OCMATHET.velocitycoord=1:length(dxdt(:,1));
OCMATHET.velocityvector=OCMATHET.velocityvector/norm(OCMATHET.velocityvector);
OCMATHET.initialpoint=depvar(:,1);
OCMATHET.differencevector=depvar2(:,end)-depvar(:,end);

% mode and path specific variables
OCMATHET.parametervalue=parametervalue(ocObj);
OCMATHET.autonomous=isautonomous(ocObj);
OCMATHET.initialtime=sol.x0;
OCMATHET.varyparameterindex=parindex;

OCMATCONT.codimension=1;

pathname=OCMATHET.datapath();
[resultfile,globalvarfile]=OCMATHET.saveintermediatefiles();
OCMATHET.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATHET.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

% free parameters
% order: switchingtimes | equilbrium | subspace parameters | free model
% parametervalues
% 

% switching times
sol.parameters=[sol.parameters switchtimes];

OCMATHET.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

% add free model parametervalues
sol.parameters=[sol.parameters 0 parametervalue(ocObj,parindex)];
OCMATHET.parametercoord=length(sol.parameters);
OCMATHET.diffperiodcoord=length(sol.parameters)-1;


function sol=generatesolstruct(ocHom,varargin)

solverInfo=solverinfo(ocHom(1));
if isfield(solverInfo,'conttype') && strcmp(solverInfo.conttype,'homoclinic')
else
    sol.parameters=[];
end
sol.x=independentvar(ocHom(1));
sol.y=dependentvar(ocHom(1));
sol.arcarg=arcargument(ocHom(1));
sol.arcinterval=arcinterval(ocHom(1));

x0=initialtime(ocHom(1));
for ii=2:2
    sol.x=[sol.x independentvar(ocHom(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocHom(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocHom(ii))];
    actarcinterval=arcinterval(ocHom(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
end
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];


function  Q=computeBase(J,stableflag,NSub)

if ~stableflag
    J=-J;
end
[VU, DU] = realeig(J);
% Select first NSub eigenvectors: unstable eigenspace
% Compute orthonormal basis for the eigenspace
VU = VU(:,1:NSub);
[Q,RU] = qr(VU);
