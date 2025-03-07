function sol=initocmat_AE_FT_IS(ocObj,ocMP,parindex,opt,varargin)
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
indifferenceorder=multiplicity(ocMP);
pathtpe=cell(1,indifferenceorder);
movinghorizon=[];
freevector=[];
hitstatevalue=[];
hitstatecoordinate=[];

for ii=1:indifferenceorder
    if isocasymptotic(ocMP(ii))
        pathtpe{ii}=pathtype(ocMP(ii));
        switch pathtpe{ii}
            case 's'
                OCMATINDIF.stableflag{ii}=1;
            case {'u','stu'}
                OCMATINDIF.stableflag{ii}=0;
        end
        OCMATINDIF.ocasymptotic(ii)=1;
    else
        OCMATINDIF.stableflag{ii}=[];
        OCMATINDIF.ocasymptotic(ii)=0;
    end
    %OCMATINDIF.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end

optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
freevectoridx=find(strcmpi(varargin,'freevector'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
simpleidx=find(strcmpi(varargin,'simple'));
movinghorizonidx=find(strcmpi(varargin,'movinghorizon'));
hitstatevalueidx=find(strcmpi(varargin,'hitstatevalue'));
hitstatecoordinateidx=find(strcmpi(varargin,'hitstatecoordinate'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
else
    targetparametervalue=[];
end
if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
else
    optimalhorizon=zeros(1,indifferenceorder);
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
else
    fixendstate=cell(1,indifferenceorder);
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
else
    fixinitstate=cell(1,indifferenceorder);
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(hitstatevalueidx)
    hitstatevalue=varargin{hitstatevalueidx+1};
end
if ~isempty(hitstatecoordinateidx)
    hitstatecoordinate=varargin{hitstatecoordinateidx+1};
end

if ~isempty(simpleidx)
    simpleflag=varargin{simpleidx+1};
else
    simpleflag=0;
end
if ~isempty(movinghorizonidx)
    movinghorizon=varargin{movinghorizonidx+1};
end
if isempty(movinghorizon)
    movinghorizon=zeros(1,indifferenceorder);
elseif numel(movinghorizon)==1
    movinghorizon=repmat(movinghorizon,1,indifferenceorder);
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if nargin==4
    opt=defaultocoptions;
end
if isempty(hitstatevalue)
    hitstatevalue=[];
    hitstatecoordinate=[];
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

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
OCMATINDIF.bcindifference=funch{5}{5};
try
    OCMATINDIF.bcoptimalhorizon=funch{5}{8};
catch
    OCMATINDIF.bcoptimalhorizon=[];
end
% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcjacobianasymptotic=funch{6}{2};
OCMATINDIF.bcjacobiantransversality=funch{6}{3};

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
OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.simple=simpleflag;
OCMATINDIF.movinghorizon=movinghorizon;
OCMATINDIF.hitstatevalue=hitstatevalue;
OCMATINDIF.hitstatecoordinate=hitstatecoordinate;

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=cell(1,indifferenceorder);
OCMATINDIF.optimalhorizon=optimalhorizon;

arcoffset=0;
limSet=cell(1,indifferenceorder);
J=cell(1,indifferenceorder);
parameters=[];
truncationtime=[];
for ii=1:indifferenceorder
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    truncationtime=[truncationtime arcintv(end)];
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    if isocasymptotic(ocMP(ii))
        limSet{ii}=limitset(ocMP(ii));
        if isequilibrium(limSet{ii})
            OCMATINDIF.limitsettype{ii}='e';
            OCMATINDIF.truncationtimecoord{ii}=[];
            limSetdepvar=dependentvar(limSet{ii});
            limSetDim=length(limSetdepvar);
            %         OCMATINDIF.truncationtimecoord{ii}=truncationtimecounter+solverInfoStruct.truncationtimecoord;
            %         truncationtimecounter=truncationtimecounter+1;
        else
            ocmaterror('Not implemented yet.')
        end
        J{ii}=linearization(limSet{ii});
        OCMATINDIF.inftimetransformation(ii)=0;
        switch OCMATINDIF.limitsettype{ii}
            case 'e'
                [asymptoticmatrix OCMATINDIF.numstable{ii} OCMATINDIF.numunstable{ii} OCMATINDIF.numcenter{ii}]=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance);
                OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
                switch pathtpe{ii}
                    case 's'
                        OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii};
                    case {'u','stu'}
                        OCMATINDIF.subspacedim{ii}=OCMATINDIF.numunstable{ii};
                    case {'sc','cs'}
                        OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii}+OCMATINDIF.numcenter{ii};
                end
                OCMATCONT.HE.equilibriumcoord{ii}=(1:length(limSetdepvar))+length(parameters);
                parameters=[parameters limSetdepvar.'];
                if simpleflag
                    OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix;
                else
                    OCMATINDIF.orthspacedim{ii}=limSetDim-OCMATINDIF.subspacedim{ii};
                    Y=zeros(OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
                    OCMATINDIF.Y{ii}=Y;
                    OCMATINDIF.Q0{ii}=computeBase(J{ii},OCMATINDIF.stableflag{ii},OCMATINDIF.subspacedim{ii});
                    OCMATCONT.HE.Ycoord{ii}=reshape(length(parameters)+(1:OCMATINDIF.orthspacedim{ii}*OCMATINDIF.subspacedim{ii}),OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
                    parameters=[parameters Y(:).'];
                    OCMATINDIF.Id{ii}=eye(OCMATINDIF.orthspacedim{ii});
                    OCMATINDIF.numY{ii}=numel(Y);
                end
        end
        OCMATINDIF.fixendstatecoord{ii}=[];
        OCMATINDIF.fixinitstatecoord{ii}=[];
        OCMATINDIF.endstate{ii}=[];
        OCMATINDIF.initstate{ii}=[];
    else
        depvar=dependentvar(ocMP(ii));
        OCMATINDIF.fixendstatecoord{ii}=fixendstate{ii};
        OCMATINDIF.fixinitstatecoord{ii}=fixinitstate{ii};
        if ~isempty(fixendstate)
            OCMATINDIF.endstate{ii}=depvar(fixendstate{ii},end);
        end
        if ~isempty(fixinitstate)
            OCMATINDIF.initstate{ii}=depvar(fixinitstate{ii},1);
        end
        OCMATINDIF.asymptoticmatrix{ii}=[];
        OCMATINDIF.orthspacedim{ii}=[];
        OCMATINDIF.Y{ii}=[];
        OCMATINDIF.Id{ii}=[];
        OCMATINDIF.numY{ii}=[];
        OCMATINDIF.inftimetransformation(ii)=0;
        limSet{ii}=[];
    end
    OCMATINDIF.numarc(ii)=arcn;
    %OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcn;
end
%OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
%OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution

counter=0;

for ii=1:indifferenceorder
    if optimalhorizon(ii)
        parameters=[parameters truncationtime(ii)];
        OCMATINDIF.optimalhorizoncoord{ii}=length(parameters);
    else
        OCMATINDIF.truncationtime(ii)=truncationtime(ii);
        OCMATINDIF.optimalhorizoncoord{ii}=[];
    end
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    if ~isempty(switchtimes{ii})
        OCMATINDIF.switchtimecoord{ii}=length(parameters)+(1:length(switchtimes{ii}));
        parameters=[parameters switchtimes{ii}];
    else
        OCMATINDIF.switchtimecoord{ii}=[];
    end
end
depvar=dependentvar(ocMP(1));

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.parameterindex=parindex;
OCMATINDIF.statecoord=statecoord(ocObj);
OCMATINDIF.costatecoord=costatecoord(ocObj);

OCMATINDIF.indifferenceorder=indifferenceorder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
%righttimeindex=cumsum(OCMATINDIF.numarc+1);
%OCMATINDIF.truncationtime=truncationtime;
OCMATINDIF.linearization=J;
OCMATINDIF.statecoordinate=1:statenum(ocObj);
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;


%OCMATINDIF.switchtimecoord=length(parameters)+(1:length(sol.parameters));
sol.parameters=[parameters sol.parameters];
OCMATINDIF.distance=[];
idx=find(movinghorizon);
if ~isempty(idx)
    OCMATINDIF.distance=zeros(1,length(movinghorizon));
    OCMATINDIF.movinghorizoncoord=zeros(1,length(movinghorizon));
    for ii=idx
        OCMATINDIF.distance(ii)=norm(OCMATINDIF.saddlepoint{ii}-ocMP(ii).y(:,end));
        OCMATINDIF.movinghorizoncoord(ii)=length(sol.parameters)+1;
        arcint=arcinterval(ocMP(ii));
        sol.parameters=[sol.parameters arcint(end)];
    end
end
OCMATINDIF.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];
OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.startvalue=depvar(OCMATINDIF.statecoord,1);

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)
global OCMATCONT
nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.y=sol.y(1:OCMATCONT.DOMAINDDATA(1).numode,:);
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
    tmpy=dependentvar(ocMultiPath(ii));
    sol.y=[sol.y tmpy(1:OCMATCONT.DOMAINDDATA(1).numode,:);];
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


