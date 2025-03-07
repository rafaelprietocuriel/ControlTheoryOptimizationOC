function sol=initocmat_EP_HET_CYC(ocObj,ocMP,parindex,fixcoordinate,opt,varargin)
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

clear global OCMATCONT OCMATHET

global OCMATCONT OCMATHET
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
hetorder=multiplicity(ocMP);
pathtpe=cell(1,hetorder);
stop=[];
userconstraint=[];
for ii=1:hetorder
    pathtpe{ii}=pathtype(ocMP(ii));
    OCMATHET.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end

targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
targetparameterindexidx=find(strcmpi(varargin,'targetparameterindex'));
connectionnumidx=find(strcmpi(varargin,'connectionnum'));
simpleidx=find(strcmpi(varargin,'simple'));
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
stopidx=find(strcmpi(varargin,'stop'));
userconstraintidx=find(strcmpi(varargin,'userconstraint'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
else
    targetparametervalue=[];
end
if ~isempty(targetparameterindexidx)
    targetparameterindex=varargin{targetparameterindexidx+1};
else
    targetparameterindex=[];
end
if ~isempty(connectionnumidx)
    connectionnum=varargin{connectionnumidx+1};
else
    connectionnum=2;
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if ~isempty(stopidx)
    stop=varargin{stopidx+1};
end
if ~isempty(userconstraintidx)
    userconstraint=varargin{userconstraintidx+1};
end

if ~isempty(simpleidx)
    simpleflag=varargin{simpleidx+1};
else
    simpleflag=0;
end
if isempty(targetparameterindex)
    targetparameterindex=parindex(end);
end
if isempty(movinghorizon)
    movinghorizon=zeros(1,hetorder);
elseif numel(movinghorizon)==1
    movinghorizon=repmat(movinghorizon,1,hetorder);
end

if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if nargin==4
    opt=defaultocoptions;
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATHET) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATHET.canonicalsystem=funch{1};
OCMATHET.canonicalsystemjacobian=funch{2}{1};
OCMATHET.canonicalsystemparameterjacobian=funch{2}{2};
OCMATHET.canonicalsystemhessian=funch{3}{1};
OCMATHET.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATHET.bcinitial=funch{5}{1};
OCMATHET.bcasymptotic=funch{5}{2};
OCMATHET.bctransversality=funch{5}{3};
OCMATHET.bcindifference=funch{5}{5};

% function for Jacobian
OCMATHET.bcjacobianinitial=funch{6}{1};
OCMATHET.bcjacobianasymptotic=funch{6}{2};
OCMATHET.bcjacobiantransversality=funch{6}{3};

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

% general function
OCMATHET.plotcontinuation=funch{11};
OCMATHET.testadmissibility=funch{12};
OCMATHET.datapath=funch{20};
OCMATHET.saveintermediatefiles=funch{21};

OCMATHET.userconstraint=userconstraint;
hybridinfo=OCMATHET.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATHET.domain(hybridinfo.arcarg(ii));
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
OCMATHET.stateconstraint=stateconstraint(ocObj);
OCMATHET.targetparametervalue=targetparametervalue;
OCMATHET.targetparameterindex=targetparameterindex;
OCMATHET.simple=simpleflag;
OCMATHET.movinghorizon=movinghorizon;
if ~isempty(stop)
    switch lower(stop)
        case 'fold'
            stop=1;
    end
end
OCMATHET.stop=stop;

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=cell(1,hetorder);
arcoffset=0;
limSet=cell(1,hetorder);
J=cell(1,hetorder);
parameters=[];
limSetdepvar0=[];
switchtimes=cell(1,hetorder);
truncationtime=zeros(1,hetorder);
for ii=1:hetorder
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    truncationtime(ii)=arcintv(end);
    arcarg=arcargument(ocMP(ii));
    OCMATHET.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    limSet{ii}=limitset(ocMP(ii));
    if isequilibrium(limSet{ii})
        OCMATHET.limitsettype{ii}='e';
        OCMATHET.truncationtimecoord{ii}=[];
        limSetdepvar=dependentvar(limSet{ii});
        limSetDim=length(limSetdepvar);
    else
        ocmaterror('Not implemented yet.')
    end
    J{ii}=linearization(limSet{ii});
    if isempty(inftimetransformation(ocMP(ii)))
        OCMATHET.inftimetransformation(ii)=0;
    else
        OCMATHET.inftimetransformation(ii)=inftimetransformation(ocMP(ii));
    end
    switch OCMATHET.limitsettype{ii}
        case 'e'
            [asymptoticmatrix OCMATHET.numstable{ii} OCMATHET.numunstable{ii} OCMATHET.numcenter{ii}]=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance);
            OCMATHET.saddlepoint{ii}=dependentvar(limSet{ii});
            switch pathtpe{ii}
                case 's'
                    OCMATHET.subspacedim{ii}=OCMATHET.numstable{ii};
                case {'u','stu'}
                    OCMATHET.subspacedim{ii}=OCMATHET.numunstable{ii};
                case {'sc','cs'}
                    OCMATHET.subspacedim{ii}=OCMATHET.numstable{ii}+OCMATHET.numcenter{ii};
            end
            if ii<=connectionnum
                OCMATHET.equilibriumcoord{ii}=(1:length(limSetdepvar))+length(parameters);
                parameters=[parameters limSetdepvar.'];
                limSetdepvar0=[limSetdepvar0 limSetdepvar];
            end
            if simpleflag
                OCMATHET.asymptoticmatrix{ii}=asymptoticmatrix;
            else
                OCMATHET.orthspacedim{ii}=limSetDim-OCMATHET.subspacedim{ii};
                Y=zeros(OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});
                OCMATHET.Y{ii}=Y;
                OCMATHET.Q0{ii}=computeBase(J{ii},OCMATHET.stableflag{ii},OCMATHET.subspacedim{ii});
                OCMATHET.Ycoord{ii}=reshape(length(parameters)+(1:OCMATHET.orthspacedim{ii}*OCMATHET.subspacedim{ii}),OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});
                parameters=[parameters Y(:).'];
                OCMATHET.Id{ii}=eye(OCMATHET.orthspacedim{ii});
                OCMATHET.numY{ii}=numel(Y);
            end
            idx=find(sqrt(sum((limSetdepvar0-repmat(limSetdepvar,1,size(limSetdepvar0,2))).^2))<5e-1);
            OCMATHET.equilibriumidx(ii)=idx;
    end
    OCMATHET.numarc(ii)=arcn;
    %OCMATHET.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATHET.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcoffset+arcn;
end
sol.parameters=parameters;
%OCMATHET.initialstateindex=cumsum(OCMATHET.initialstateindex);
%OCMATHET.initialstateindex=[0 OCMATHET.initialstateindex(end-1)]+1;
OCMATHET.solutionindex=zeros(1,sum(OCMATHET.numarc));% relate arcindex to indifference solution
counter=0;

for ii=1:hetorder
    counter_start=counter+1;
    counter=counter+OCMATHET.numarc(ii);
    OCMATHET.solutionindex(counter_start:counter)=ii;
    if ii==1
        if ~isempty(switchtimes{ii})
            OCMATHET.switchtimecoord{ii}=(1:length(switchtimes{ii}))+length(parameters);
            addval=OCMATHET.switchtimecoord{ii}(end);
            sol.parameters=[sol.parameters switchtimes{ii}];
        else
            OCMATHET.switchtimecoord{ii}=[];
            addval=length(parameters);
        end
    else
        if ~isempty(switchtimes{ii})
            OCMATHET.switchtimecoord{ii}=addval+(1:length(switchtimes{ii}));
            addval=OCMATHET.switchtimecoord{ii}(end);
            sol.parameters=[sol.parameters switchtimes{ii}];
        else
            OCMATHET.switchtimecoord{ii}=[];
        end
    end
end
%depvar=dependentvar(ocMP(1));

OCMATHET.cumsumnumarc=cumsum(OCMATHET.numarc);
OCMATHET.initcoord=[1 OCMATHET.cumsumnumarc(1:end-1)+1];
OCMATHET.parameterindex=parindex;
OCMATHET.initialstatedifference=zeros(length(ocMP(1).y(:,1)),connectionnum);
OCMATHET.fixvalue=zeros(length(fixcoordinate),connectionnum);
for ii=1:connectionnum
    OCMATHET.initialstatedifference(:,ii)=ocMP(2*ii).y(:,1)-ocMP(2*ii-1).y(:,1);
    OCMATHET.fixvalue(:,ii)=ocMP(2*ii).y(fixcoordinate,1);
end
OCMATHET.hetorder=hetorder;
OCMATHET.connectionnum=connectionnum;
OCMATHET.parametervalue=parametervalue(ocObj);
OCMATHET.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATHET.numarc(1:end-1)+1)];
%righttimeindex=cumsum(OCMATHET.numarc+1);
OCMATHET.truncationtime=truncationtime;
OCMATHET.linearization=J;
OCMATHET.statecoordinate=1:statenum(ocObj);
OCMATHET.fixcoordinate=fixcoordinate;
OCMATHET.pathtype=pathtpe;
pathname=OCMATHET.datapath();
[resultfile,globalvarfile]=OCMATHET.saveintermediatefiles();
OCMATHET.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATHET.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATHET.objectivevaluecalc=0;
OCMATHET.autonomous=1;
OCMATHET.distance=[];
idx=find(movinghorizon);
if ~isempty(idx)
    OCMATHET.distance=zeros(1,length(movinghorizon));
    OCMATHET.movinghorizoncoord=zeros(1,length(movinghorizon));
    for ii=idx
        OCMATHET.distance(ii)=norm(OCMATHET.saddlepoint{ii}-ocMP(ii).y(:,end));
        OCMATHET.movinghorizoncoord(ii)=length(sol.parameters)+1;
        arcint=arcinterval(ocMP(ii));
        sol.parameters=[sol.parameters arcint(end)];
    end
end
OCMATHET.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

sol.parameters=[sol.parameters OCMATHET.parametervalue(parindex)];

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
%sol.parameters=parameters(ocMultiPath(1));
%numcontpar=length(continuationparameter(ocMultiPath(1)));
%sol.parameters(end-numcontpar+1:end)=[];
%if isempty(sol.parameters)
%    sol.parameters=sol.arcinterval(2:end-1);
%end
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
    %    freepar=parameters(ocMultiPath(ii));
    %    if isempty(freepar)
    %        freepar=actarcinterval(2:end-1);
    %        numcontpar=0;
    %    else
    %        numcontpar=length(continuationparameter(ocMultiPath(ii)));
    
    %    end
    %    freepar(end-numcontpar+1:end)=[];
    %    sol.parameters=[sol.parameters freepar];
end
sol.parameters=[];
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


