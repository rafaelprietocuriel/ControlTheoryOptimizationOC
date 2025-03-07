function sol=initocmat_AE_IS_P(ocObj,ocMP,parindex,opt,varargin)
%
% INITOCMAT_AE_IS initialization for the continuation of an indifference
% threshold
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE) the
% continuation of an indifference threshold in the state space is
% initialized.
% OCOBJ          ... corresponding optimal control model
% OCMP           ... two cell array of ocasymptotics, for the different
%                    solution paths or an instance of an ocmultipath
%                    object.
% TARGETCOORDINATE ... the continuation is done along the n-1 coordinates
%                   (n number of states)
% TARGETVALUE    ... The value of the target vector for the n-1
%                   coordinates.
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and OCMATINDIF
% are initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

clear global OCMATCONT OCMATINDIF

global OCMATCONT OCMATINDIF
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
stopcriterion=[];
freevector=[];
freeendtime=[];
targetparametervalue=[];
targetparameterindex=[];
for ii=1:indifforder
    pathtpe{ii}=pathtype(ocMP(ii));
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end

stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
targetparameterindexidx=find(strcmpi(varargin,'targetparameterindex'));
freevectoridx=find(strcmpi(varargin,'freevector'));
excludecoordinate4epidx=find(strcmpi(varargin,'excludecoordinate4ep'));
jacobianidx=find(strcmpi(varargin,'jacobian'));
arcid4epidx=find(strcmpi(varargin,'arcid4ep'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(jacobianidx)
    jacobian=varargin{jacobianidx+1};
else
    jacobian=cell(1,indifforder);
end
if ~isempty(excludecoordinate4epidx)
    excludecoordinate4ep=varargin{excludecoordinate4epidx+1};
else
    excludecoordinate4ep=cell(1,indifforder);
end
if ~isempty(arcid4epidx)
    arcid4ep=varargin{arcid4epidx+1};
else
    arcid4ep=cell(1,indifforder);
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};  
end
if ~isempty(targetparameterindexidx)
    targetparameterindex=varargin{targetparameterindexidx+1};  
    if ischar(targetparameterindex)
        targetparameterindex=parameterindex(ocObj,targetparameterindex);
    end
end
if isempty(stopcriterion)
    stopcriterion=0;
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder); % stop if limitpoint occurs
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(targetparameterindex) && length(parindex)==1
    targetparameterindex=parindex;
end
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end
OCMATINDIF.stopcriterion=stopcriterion;

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     ocMP(ii).solver='';
%     ocMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4IndifferenceSolutionContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATINDIF.canonicalsystem=funch{1};
OCMATINDIF.canonicalsystemjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
OCMATINDIF.canonicalsystemhessian=funch{3}{1};
OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.equilibrium=funch{5}{4};
OCMATINDIF.bcindifference=funch{5}{5};

% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcjacobianasymptotic=funch{6}{2};
OCMATINDIF.bcjacobiantransversality=funch{6}{3};
OCMATINDIF.bcjacobianindifference=funch{6}{4};

% function describing the hybrid structure of the problem
OCMATINDIF.hybridinfo=funch{7}{1};
OCMATINDIF.domain=funch{7}{2};
OCMATINDIF.guard=funch{7}{3};
OCMATINDIF.reset=funch{7}{4};
OCMATINDIF.switchtime=funch{7}{5};
OCMATINDIF.jacobianguard=funch{7}{7};
OCMATINDIF.jacobianreset=funch{7}{8};
OCMATINDIF.domaindiscretization=funch{7}{9};
OCMATINDIF.timesettransformation=funch{7}{10};

% general function
OCMATINDIF.plotcontinuation=funch{11};
OCMATINDIF.testadmissibility=funch{12};
OCMATINDIF.datapath=funch{20};
OCMATINDIF.saveintermediatefiles=funch{21};

hybridinfo=OCMATINDIF.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATINDIF.domain(hybridinfo.arcarg(ii));
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
% test if pure state constraints are defined
OCMATINDIF.stateconstraint=stateconstraint(ocObj);
OCMATINDIF.freeendtime=freeendtime;
OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.targetparameterindex=targetparameterindex;

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
limSet=cell(1,indifforder);
J=cell(1,indifforder);
parameters=[];
truncationtime=[];
for ii=1:indifforder
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    truncationtime=[truncationtime arcintv(end)];
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    limSet{ii}=limitset(ocMP(ii));
    OCMATINDIF.limitsettype{ii}='e';
    OCMATINDIF.truncationtimecoord{ii}=[];
    limSetdepvar=dependentvar(limSet{ii});
    if ~isempty(excludecoordinate4ep{ii})
        limSetdepvar(excludecoordinate4ep{ii})=[];
    end
    limSetDim=length(limSetdepvar);
    if ~isempty(jacobian{ii})
        J{ii}=jacobian{ii};
    else
        J{ii}=linearization(limSet{ii});
        J{ii}(:,excludecoordinate4ep{ii})=[];
        J{ii}(excludecoordinate4ep{ii},:)=[];
    end
    if ~isempty(arcid4ep{ii})
        OCMATINDIF.arcid4ep(ii)=arcid4ep{ii};
    else
        OCMATINDIF.arcid4ep(ii)=arcargument(limSet{ii});
    end
    
    if isempty(inftimetransformation(ocMP(ii)))
        OCMATINDIF.inftimetransformation(ii)=0;
    else
        OCMATINDIF.inftimetransformation(ii)=inftimetransformation(ocMP(ii));
    end
    [dum OCMATINDIF.numstable{ii} OCMATINDIF.numunstable{ii} OCMATINDIF.numcenter{ii}]=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance);
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
    switch pathtpe{ii}
        case {'sts','s'}
            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii};
        case 'u'
            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numunstable{ii};
        case 'wu'
            OCMATINDIF.subspacedim{ii}=1;
        case {'sc','cs'}
            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii}+OCMATINDIF.numcenter{ii};
    end
    OCMATINDIF.orthspacedim{ii}=limSetDim-OCMATINDIF.subspacedim{ii};
    Y=zeros(OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
    OCMATINDIF.Y{ii}=Y;
    OCMATINDIF.Q0{ii}=computeBase(J{ii},1,OCMATINDIF.subspacedim{ii});
    OCMATCONT.HE.equilibriumcoord{ii}=(1:length(limSetdepvar))+length(parameters);
    parameters=[parameters limSetdepvar.'];
    OCMATCONT.HE.Ycoord{ii}=reshape(length(parameters)+(1:OCMATINDIF.orthspacedim{ii}*OCMATINDIF.subspacedim{ii}),OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
    parameters=[parameters Y(:).'];
    OCMATINDIF.Id{ii}=eye(OCMATINDIF.orthspacedim{ii});
    OCMATINDIF.numY{ii}=numel(Y);
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    if freeendtime(ii)>0
        depvar=dependentvar(ocMP(ii));
        OCMATINDIF.distance{ii}=norm(OCMATINDIF.saddlepoint{ii}-depvar(:,end));
    else
        OCMATINDIF.distance{ii}=[];
    end
    arcoffset=arcoffset+arcn;
end
OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution
sol.parameters=[parameters switchtimes{:}];
counter=0;
paroffset=length(parameters);
for ii=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    if ~isempty(switchtimes{ii})
        OCMATINDIF.switchtimecoord{ii}=paroffset+(1:length(switchtimes{ii}));
    else
        OCMATINDIF.switchtimecoord{ii}=[];
    end
    if ~isempty(OCMATINDIF.switchtimecoord{ii})
        paroffset=OCMATINDIF.switchtimecoord{ii}(end);
    end
end
righttimeindex=cumsum(OCMATINDIF.numarc+1);
OCMATINDIF.truncationtimecoord=[];
ctr=0;
for ii=1:indifforder
    if freeendtime(ii)
        ctr=ctr+1;
        OCMATINDIF.truncationtimecoord{ii}=paroffset+ctr;
        sol.parameters=[sol.parameters sol.arcinterval(righttimeindex(ii))];
        OCMATINDIF.truncationtime{ii}=[];
    else
        OCMATINDIF.truncationtime{ii}=sol.arcinterval(righttimeindex(ii));
    end
end
OCMATINDIF.statecoordinate=1:statenum(ocObj);

depvar=dependentvar(ocMP(1));
OCMATINDIF.startvalue=depvar(OCMATINDIF.statecoordinate,1);
OCMATINDIF.excludecoordinate4ep=excludecoordinate4ep;

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.parameterindex=parindex;

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
righttimeindex=cumsum(OCMATINDIF.numarc+1);
OCMATINDIF.truncationtime=sol.arcinterval(righttimeindex);
OCMATINDIF.linearization=J;
OCMATINDIF.freevector=freevector;
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;

%sol.parameters=parameters;
OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.parametervaluecoord=length(sol.parameters)+(1:length(parindex));
if length(parindex)>1
    OCMATINDIF.targetparametervaluecoord=OCMATINDIF.parametervaluecoord(parindex==targetparameterindex);
else
    OCMATINDIF.targetparametervaluecoord=OCMATINDIF.parametervaluecoord;
end
sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
sol.parameters=parameters(ocMultiPath(1));
%numcontpar=length(continuationparameter(ocMultiPath(1)));
% sol.parameters(end-numcontpar+1:end)=[];
% if isempty(sol.parameters)
%     sol.parameters=sol.arcinterval(2:end-1);
% end
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
%     freepar=parameters(ocMultiPath(ii));
%     if isempty(freepar)
%         freepar=actarcinterval(2:end-1);
%         numcontpar=0;
%     else
%         numcontpar=length(continuationparameter(ocMultiPath(ii)));
% 
%     end
%     freepar(end-numcontpar+1:end)=[];
%     sol.parameters=[sol.parameters freepar];
end
sol.parameters=[];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
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


