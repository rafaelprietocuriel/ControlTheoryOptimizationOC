function sol=initocmat_AE_Hom(ocObj,ocMP,parindex,opt,varargin)
% INITOCMAT_LC_H_P initialization for the continuation of a limit cycle
% with respect to a parameter, starting from a Hopf bifurcation
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX)
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATHOM
global OCMATCONT OCMATHOM
sol=[];

if nargin==3
    opt=defaultocoptions;
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc limitcycle is empty.')
    return
end
if nargin==4
    opt=[];
end
if isempty(opt)
    opt=defaultocoptions;
end
ocMP=ocmultipath(ocMP);
for ii=1:2
    pathtpe{ii}=pathtype(ocMP(ii));
end
if ~sum(strcmp(pathtpe,'s'))==1  ||  ~sum(strcmp(pathtpe,'u'))==1
    ocmaterror('Input multipath has to consist of a stable and unstable path.')
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

OCMATHOM.canonicalsystem=funch{1};
OCMATHOM.canonicalsystemjacobian=funch{2}{1};
OCMATHOM.canonicalsystemparameterjacobian=funch{2}{2};
OCMATHOM.canonicalsystemhessian=funch{3}{1};
OCMATHOM.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions

% function for Jacobian

% function describing the hybrid structure of the problem
OCMATHOM.hybridinfo=funch{7}{1};
OCMATHOM.domain=funch{7}{2};
OCMATHOM.guard=funch{7}{3};
OCMATHOM.reset=funch{7}{4};
OCMATHOM.switchtime=funch{7}{5};
OCMATHOM.jacobianguard=funch{7}{7};
OCMATHOM.jacobianreset=funch{7}{8};
OCMATHOM.domaindiscretization=funch{7}{9};
OCMATHOM.timesettransformation=funch{7}{10};
if ~isautonomous(ocObj)
    OCMATHOM.canonicalsystemderivativetime=funch{2}{3};
end

% function for the boundary conditions
OCMATHOM.bcasymptotic=funch{5}{2};
% general function
OCMATHOM.plotcontinuation=funch{11};
OCMATHOM.testadmissibility=funch{12};
OCMATHOM.datapath=funch{20};
OCMATHOM.saveintermediatefiles=funch{21};

hybridinfo=OCMATHOM.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATHOM.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocMP);
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% limit set
limSet=limitset(ocMP(1));
J=linearization(limSet);
limSetdepvar=dependentvar(limSet);
limSetDim=length(limSetdepvar);

% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
switchtimes=[];
for ii=1:2
    arcn=arcnum(ocMP(ii));
    %solverInfoStruct=solverinfo(ocMP(ii));
    OCMATHOM.switchtimecoord{ii}=length(switchtimes)+(1:arcn-1);
    switchtimes=[switchtimes ocMP(ii).arcinterval(2:end-1)];
    arcarg=arcargument(ocMP(ii));

    OCMATHOM.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATHOM.inftimetransformation(ii)=0;
    [OCMATHOM.asymptoticmatrix{ii},OCMATHOM.numstable,OCMATHOM.numunstable,OCMATHOM.numcenter]=asymptoticbc(J,pathtpe{ii},'c',ZeroDeviationTolerance);
    OCMATHOM.numarc(ii)=arcn;
    OCMATHOM.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATHOM.arccoord{ii}=[1:arcn]+arcoffset;
    arcoffset=arcn;
end
OCMATHOM.saddlepoint=dependentvar(limSet);
OCMATHOM.pathtype=pathtpe;
OCMATHOM.initialstateindex=cumsum(OCMATHOM.initialstateindex);
OCMATHOM.initialstateindex=[0 OCMATHOM.initialstateindex(end-1)]+1;
OCMATHOM.solutionindex=zeros(1,sum(OCMATHOM.numarc));% relate arcindex to indifference solution
counter=0;
for ii=1:2
    counter_start=counter+1;
    counter=counter+OCMATHOM.numarc(ii);
    OCMATHOM.solutionindex(counter_start:counter)=ii;
end

OCMATHOM.cumsumnumarc=cumsum(OCMATHOM.numarc);
OCMATHOM.initcoord=[1 OCMATHOM.cumsumnumarc(1:end-1)+1];

OCMATHOM.parametervalue=parametervalue(ocObj);
OCMATHOM.initialtime=sol.x0;

righttimeindex=cumsum(OCMATHOM.numarc+1);
OCMATHOM.truncationtime=sol.arcinterval(righttimeindex);

% phase condition
dxdt=canonicalsystem(ocObj,ocMP(1),[],1);
depvar=dependentvar(ocMP(1));
OCMATHOM.velocityvector=dxdt(:,1);
OCMATHOM.velocitycoord=1:length(dxdt(:,1));
OCMATHOM.velocityvector=OCMATHOM.velocityvector/norm(OCMATHOM.velocityvector);
OCMATHOM.initialpoint=depvar(:,1);

% subspace continuation
for ii=1:2
    stableflag=~isempty(strfind(pathtpe{ii},'s'));
    switch pathtpe{ii}
        case 's'
            OCMATHOM.subspacedim{ii}=OCMATHOM.numstable;
        case 'u'
            OCMATHOM.subspacedim{ii}=OCMATHOM.numunstable;
    end
    OCMATHOM.orthspacedim{ii}=limSetDim-OCMATHOM.subspacedim{ii};
    Y=zeros(OCMATHOM.orthspacedim{ii},OCMATHOM.subspacedim{ii});
    OCMATHOM.Y{ii}=Y;
    OCMATHOM.Q0{ii}=computeBase(J,stableflag,OCMATHOM.subspacedim{ii});
    OCMATHOM.Id{ii}=eye(OCMATHOM.orthspacedim{ii});
    OCMATHOM.numY{ii}=numel(Y);
end
% mode and path specific variables
OCMATHOM.parametervalue=parametervalue(ocObj);
OCMATHOM.autonomous=isautonomous(ocObj);
OCMATHOM.initialtime=sol.x0;
OCMATHOM.varyparameterindex=parindex;

OCMATCONT.codimension=1;

pathname=OCMATHOM.datapath();
[resultfile,globalvarfile]=OCMATHOM.saveintermediatefiles();
OCMATHOM.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATHOM.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

% free parameters
% order: switchingtimes | equilbrium | subspace parameters | free model
% parametervalues
% 

% switching times
sol.parameters=[sol.parameters switchtimes];
% saddle point
sol.parameters=[sol.parameters OCMATHOM.saddlepoint(:).'];
OCMATHOM.equilibriumcoord=length(sol.parameters)-length(OCMATHOM.saddlepoint)+1:length(sol.parameters);

% subspace parameters
for ii=1:2
    Y=OCMATHOM.Y{ii};
    OCMATCONT.HE.Ycoord{ii}=reshape(length(sol.parameters)+(1:OCMATHOM.orthspacedim{ii}*OCMATHOM.subspacedim{ii}),OCMATHOM.orthspacedim{ii},OCMATHOM.subspacedim{ii});
    sol.parameters=[sol.parameters Y(:).'];
end
OCMATHOM.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

% add free model parametervalues
sol.parameters=[sol.parameters parametervalue(ocObj,parindex)];


% total 
OCMATHOM.parametercoord=1:length(sol.parameters);

function sol=generatesolstruct(ocMP,varargin)

solverInfo=solverinfo(ocMP(1));
if isfield(solverInfo,'conttype') && strcmp(solverInfo.conttype,'homoclinic')
else
    sol.parameters=[];
end
sol.x=independentvar(ocMP(1));
sol.y=dependentvar(ocMP(1));
sol.arcarg=arcargument(ocMP(1));
sol.arcinterval=arcinterval(ocMP(1));

x0=initialtime(ocMP(1));
for ii=2:2
    sol.x=[sol.x independentvar(ocMP(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMP(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMP(ii))];
    actarcinterval=arcinterval(ocMP(ii));
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
