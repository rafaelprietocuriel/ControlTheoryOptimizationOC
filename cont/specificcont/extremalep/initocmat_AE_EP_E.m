function sol=initocmat_AE_EP_E(ocObj,ocEP,contcoordinate,targetvalue,initStruct,opt,varargin)
%
% INITOCMAT_AE_AE_E initialization for exact asymptotic extremal
% calculation, i.e. the dynamics corresponding to the equilibrium is linear
%
% SOL=INITOCMAT_AE_AE_E(OCOBJ,INITSTRUCT,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_AE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
pathtpe='';
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(initStruct)
    ocmatmsg('init structure is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
pathtypeidx=find(strcmpi(varargin,'pathtype'));
if ~isempty(pathtypeidx)
    pathtpe=varargin{pathtypeidx+1};
end
if isempty(pathtpe)
    pathtpe='s';
end
targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcexactasymptotic=funch{5}{2};

% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianexactasymptotic=funch{6}{2};

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

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

% test if pure state constraints are defined
OCMATAE.stateconstraint=stateconstraint(ocObj);
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
J=linearization(ocEP);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=0;
OCMATAE.inftimetransformation=0;
OCMATAE.endarcarg=arcargument(ocEP);

OCMATAE.linearization=J;
OCMATAE.saddlepoint=dependentvar(ocEP);
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=initStruct.y(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;


sol=generatesolstruct(ocObj,initStruct,opt);
sol.parameters=[sol.parameters 0];

function sol=generatesolstruct(ocObj,initStruct,opt)
global OCMATAE
sol=[];
equationsolver=str2func(opt.GENERAL.EquationSolver);

numdepvar=size(initStruct.y,1);
OCMATAE.exactendpointcoord=1:numdepvar;
OCMATAE.switchtimecoord=numdepvar+1;

% test for state constraints
sccounter=0;
if OCMATAE.stateconstraint
    arcarg4sc=arcargwithactivestateconstraint(ocObj);
    if ismember(OCMATAE.endarcarg,arcarg4sc)
        sccounter=1;
        OCMATAE.jumpcostateindex=2;
    end
    if sccounter
        OCMATAE.jumpcostatecoord=OCMATAE.switchtimecoord(end)+1;
        if isfield(initStruct,'jump')
            jump=initStruct.jump;
        else
            jump=0;
        end
    else
        OCMATAE.jumpcostatecoord=[];
        OCMATAE.jumpcostateindex=[];
        jump=[];
    end
end

switchtime=initStruct.x0; % the switchtime is superfluous (because it is assumed that it is zero), but kept to easily extend the equations to the more general  case with switchtime unequal zero. 
[totaldepvar,fval]=equationsolver(@switchpointcondition,[initStruct.y(:);switchtime;jump],opt.EQ,[initStruct.arcarg;OCMATAE.endarcarg],OCMATAE.parametervalue);
if norm(fval)>opt.GENERAL.ZeroDeviationTolerance
    ocmatmsg('Detected solution may be inaccurate.')
end
if OCMATAE.stateconstraint
    jump=totaldepvar(end);
    totaldepvar(end)=[];
else
    jump=[];
end
switchtime=totaldepvar(end);
totaldepvar(end)=[];
totaldepvar=reshape(totaldepvar,[],2);
sol.y=totaldepvar(:,ones(1,opt.GENERAL.TrivialArcMeshNum));
sol.x=linspace(0,1,opt.GENERAL.TrivialArcMeshNum);
sol.x0=initStruct.x0;
sol.arcarg=initStruct.arcarg;
sol.arcposition=[1;opt.GENERAL.TrivialArcMeshNum];
sol.arcinterval=[0 switchtime];
sol.parameters=[totaldepvar(:,2).' switchtime jump];
sol.solver='';
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];

OCMATAE.arcarg=[initStruct.arcarg OCMATAE.endarcarg];
OCMATAE.numarc=2;
OCMATAE.edge=[initStruct.arcarg;OCMATAE.endarcarg];

function res=switchpointcondition(totaldepvar,edge,modelpar,varargin)
global OCMATAE
if OCMATAE.stateconstraint
    jump=totaldepvar(end);
    jump=[0 jump];
    totaldepvar(end)=[];
end
switchtime=totaldepvar(end);
totaldepvar(end)=[];
totaldepvar=reshape(totaldepvar,[],2);
depvara=[totaldepvar(:,1) totaldepvar(:,2)];
depvarb=[totaldepvar(:,1) totaldepvar(:,2)];
resinit=OCMATAE.bcinitial(depvara,OCMATAE.targetcoordinate,OCMATAE.startvalue,modelpar,edge(1));
resasym=OCMATAE.bcexactasymptotic(depvarb,OCMATAE.asymptoticmatrix,OCMATAE.saddlepoint);
if OCMATAE.stateconstraint
    resconnec=[OCMATAE.reset(depvara,depvarb,modelpar,switchtime,jump,edge,edge,1); ...
        OCMATAE.guard(depvara,depvarb,modelpar,switchtime,jump,edge,edge,1)];
else
    resconnec=OCMATAE.reset(depvara,depvarb,modelpar,0,edge,edge,1);
end
res=[resinit;resconnec;resasym];