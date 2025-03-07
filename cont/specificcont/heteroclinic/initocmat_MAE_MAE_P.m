function sol=initocmat_MAE_MAE_P(ocObj,ocMP,parindex,initialcoordinate,opt,varargin)
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
% During the initialization two global variables OCMATCONT and OCMATAE
% are initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

clear global OCMATCONT OCMATAE

global OCMATCONT OCMATAE
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
multorder=multiplicity(ocMP);
pathtpe=cell(1,multorder);
for ii=1:multorder
    pathtpe{ii}=pathtype(ocMP(ii));
    OCMATAE.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
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
simpleidx=find(strcmpi(varargin,'simple'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
else
    targetparametervalue=[];
end

if ~isempty(simpleidx)
    simpleflag=varargin{simpleidx+1};
else
    simpleflag=0;
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

% initialize global variable (OCMATAE) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.bcindifference=funch{5}{5};

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
OCMATAE.targetparametervalue=targetparametervalue;
OCMATAE.simple=simpleflag;

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=cell(1,multorder);
arcoffset=0;
parameters=[];
truncationtime=[];
for ii=1:multorder
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    truncationtime=[truncationtime arcintv(end)];
    arcarg=arcargument(ocMP(ii));
    OCMATAE.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    limSet=limitset(ocMP(ii));
    if isequilibrium(limSet)
        OCMATAE.limitsettype='e';
        OCMATAE.truncationtimecoord{ii}=[];
        limSetdepvar=dependentvar(limSet);
        limSetDim=length(limSetdepvar);
        %         OCMATAE.truncationtimecoord{ii}=truncationtimecounter+solverInfoStruct.truncationtimecoord;
        %         truncationtimecounter=truncationtimecounter+1;
    else
        ocmaterror('Not implemented yet.')
    end
    J=linearization(limSet);
    if isempty(inftimetransformation(ocMP(ii)))
        OCMATAE.inftimetransformation(ii)=0;
    else
        OCMATAE.inftimetransformation(ii)=inftimetransformation(ocMP(ii));
    end
    switch OCMATAE.limitsettype
        case 'e'
            [asymptoticmatrix OCMATAE.numstable{ii} OCMATAE.numunstable{ii} OCMATAE.numcenter{ii}]=asymptoticbc(J,pathtpe{ii},'c',ZeroDeviationTolerance);
            OCMATAE.saddlepoint=dependentvar(limSet);
            switch pathtpe{ii}
                case 's'
                    OCMATAE.subspacedim{ii}=OCMATAE.numstable{ii};
                case 'u'
                    OCMATAE.subspacedim{ii}=OCMATAE.numunstable{ii};
                case {'sc','cs'}
                    OCMATAE.subspacedim{ii}=OCMATAE.numstable{ii}+OCMATAE.numcenter{ii};
            end
            if ii==1
                OCMATCONT.HE.equilibriumcoord=(1:length(limSetdepvar))+length(parameters);
                parameters=[parameters limSetdepvar.'];
            end
            if simpleflag
                OCMATAE.asymptoticmatrix{ii}=asymptoticmatrix;
            else
                OCMATAE.orthspacedim{ii}=limSetDim-OCMATAE.subspacedim{ii};
                Y=zeros(OCMATAE.orthspacedim{ii},OCMATAE.subspacedim{ii});
                OCMATAE.Y{ii}=Y;
                OCMATAE.Q0{ii}=computeBase(J,OCMATAE.stableflag{ii},OCMATAE.subspacedim{ii});
                OCMATCONT.HE.Ycoord{ii}=reshape(length(parameters)+(1:OCMATAE.orthspacedim{ii}*OCMATAE.subspacedim{ii}),OCMATAE.orthspacedim{ii},OCMATAE.subspacedim{ii});
                parameters=[parameters Y(:).'];
                OCMATAE.Id{ii}=eye(OCMATAE.orthspacedim{ii});
                OCMATAE.numY{ii}=numel(Y);
            end
    end
    OCMATAE.numarc(ii)=arcn;
    %OCMATAE.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATAE.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcn;
end
%OCMATAE.initialstateindex=cumsum(OCMATAE.initialstateindex);
%OCMATAE.initialstateindex=[0 OCMATAE.initialstateindex(end-1)]+1;
OCMATAE.solutionindex=zeros(1,sum(OCMATAE.numarc));% relate arcindex to indifference solution
counter=0;

for ii=1:multorder
    counter_start=counter+1;
    counter=counter+OCMATAE.numarc(ii);
    OCMATAE.solutionindex(counter_start:counter)=ii;
    if ii==1
        if ~isempty(switchtimes{ii})
            OCMATAE.switchtimecoord{ii}=(1:length(switchtimes{ii}))+length(parameters);
            addval=OCMATAE.switchtimecoord{ii}(end);
        else
            OCMATAE.switchtimecoord{ii}=[];
            addval=length(parameters);
        end
    else
        OCMATAE.switchtimecoord{ii}=addval+(1:length(switchtimes{ii}));
    end
    OCMATAE.initialstateindex(ii)=numel(ocMP(ii).x);
end
depvar=dependentvar(ocMP(1));

OCMATAE.initialstateindex=cumsum(OCMATAE.initialstateindex);
OCMATAE.initialstateindex=[0 OCMATAE.initialstateindex(end-1)]+1;
OCMATAE.cumsumnumarc=cumsum(OCMATAE.numarc);
OCMATAE.initcoord=[1 OCMATAE.cumsumnumarc(1:end-1)+1];
OCMATAE.parameterindex=parindex;

OCMATAE.multorder=multorder;
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATAE.numarc(1:end-1)+1)];
%righttimeindex=cumsum(OCMATAE.numarc+1);
OCMATAE.truncationtime=truncationtime;
OCMATAE.linearization=J;
OCMATAE.statecoordinate=1:statenum(ocObj);
OCMATAE.initialcoordinate=initialcoordinate;
OCMATAE.initialstate=depvar(initialcoordinate,1);
OCMATAE.pathtype=pathtpe;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATAE.objectivevaluecalc=0;
OCMATAE.autonomous=1;


%OCMATAE.switchtimecoord=length(parameters)+(1:length(sol.parameters));
sol.parameters=[parameters [switchtimes{:}] sol.parameters];
OCMATAE.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

sol.parameters=[sol.parameters OCMATAE.parametervalue(parindex)];
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


