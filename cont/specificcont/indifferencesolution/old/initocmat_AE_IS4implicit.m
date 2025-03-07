function sol=initocmat_AE_IS4implicit(ocObj,ocMP,varargin)
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

global OCMATCONT OCMATINDIF
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
targetvalue=[];
stopcriterion=[];
freeendtime=[];
fixdistance=[];
fixdistancecoordinate=[];
fixendcoordinate=[];
targetcoordinate=[];
freevector=[];
targetvectorcoordinate=[];
asymptoticmatrix=[];
option=[];
divergingcoordinate=[];
exogenousfunction=[];
exogenousinitialstates=[];
userfunction=[];
hitfunction=[];

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
if nargin==2
    option=defaultocoptions;
end

for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end

if isempty(stopcriterion)
    stopcriterion=1; % stop if limitpoint occurs
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder);
end
if isempty(asymptoticmatrix)
    asymptoticmatrix=cell(1,indifforder);
end
OCMATINDIF.stopcriterion=stopcriterion;
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end
for ii=1:indifforder
    limSet{ii}=limitset(ocMP(ii));
%     if isperiodic(limSet{ii}) && strcmp(pathtpe{ii},'s')
%         freeendtime(ii)=-1;
%     end
end
if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=zeros(1,indifforder);
end
if ~isempty(divergingcoordinate)
    if ~iscell(divergingcoordinate)
        divergingcoordinate=repmat({divergingcoordinate},1,indifforder);
    end
end
if ~isempty(fixendcoordinate)
    if ~iscell(fixendcoordinate)
        fixendcoordinate=repmat({fixendcoordinate},1,indifforder);
    end
else
    fixendcoordinate=cell(1,indifforder);
end
OCMATINDIF.userfunction=~isempty(userfunction);
OCMATINDIF.hitfunction=~isempty(hitfunction);

OCMATINDIF.exogenousfunction=exogenousfunction;
if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end

ZeroDeviationTolerance=zeros(1,indifforder);
for ii=1:indifforder
    if iscell(option)
        ZeroDeviationTolerance(ii)=getocoptions(option{ii},'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
    else
        ZeroDeviationTolerance(ii)=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
    end
end
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

if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamics=funch{4}{1};
    OCMATINDIF.exogenousjacobian=funch{4}{2};
    OCMATINDIF.exogenousparameterjacobian=funch{4}{3};
    try
        OCMATINDIF.exogenousderivativetime=funch{4}{4};
    catch
        OCMATINDIF.exogenousderivativetime=[];
    end
    OCMATINDIF.exogenousinitialstatesfunc=funch{4}{11};
end

% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
if ~isempty(divergingcoordinate)
    OCMATINDIF.bcinf=funch{5}{4};
end
OCMATINDIF.bcindifference=funch{5}{5};
OCMATINDIF.algebraicequation=funch{5}{7};

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

if OCMATINDIF.hitfunction
    OCMATINDIF.targetfunction=funch{10}{1};
end
if OCMATINDIF.userfunction
    OCMATINDIF.userfunction=funch{10}{2};
end

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
OCMATINDIF.divergingcoordinate=divergingcoordinate;
OCMATINDIF.fixendcoordinate=fixendcoordinate;
OCMATINDIF.fixdistance=fixdistance;

OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.targetvectorcoordinate=targetvectorcoordinate;

OCMATINDIF.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end


scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;

if ~isempty(divergingcoordinate)
    for ii=1:indifforder
        OCMATINDIF.convergingcoordinate{ii}=setdiff(OCMATINDIF.statecostatecoord,OCMATINDIF.divergingcoordinate{ii});
    end
else
     OCMATINDIF.divergingcoordinate=cell(indifforder,1);
end

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=[];
limSet=cell(1,indifforder);
J=cell(1,indifforder);
switchtimes=cell(1,indifforder);
arcoffset=0;
truncationtimecounter=0;
numberofodes=[];
for ii=1:indifforder
    numberofodes=[numberofodes odenumber(ocMP(ii))];

    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);

    limSet{ii}=limitset(ocMP(ii));
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    if isequilibrium(limSet{ii})
        OCMATINDIF.limitsettype{ii}='e';
    elseif isperiodic(limSet{ii})
        OCMATINDIF.limitsettype{ii}='l';
        truncationtimecounter=truncationtimecounter+1;
    else
        ocmaterror('Limit set type unknown.')
    end
    tmp=jacobian(limSet{ii});
    J{ii}=tmp{1};
    switch OCMATINDIF.limitsettype{ii}
        case 'e'
            if isempty(asymptoticmatrix{ii})
                Jred=J{ii};
                if ~isempty( OCMATINDIF.divergingcoordinate{ii})
                    Jred(OCMATINDIF.divergingcoordinate{ii},:)=[];
                    Jred(:,OCMATINDIF.divergingcoordinate{ii})=[];
                end
                OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(Jred,pathtpe{ii},'c',ZeroDeviationTolerance(ii));
            else
                OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix{ii};
            end
            tmp=dependentvar(limSet{ii});
            OCMATINDIF.saddlepoint{ii}=tmp{1};
            if ~isempty(OCMATINDIF.fixendcoordinate{ii})
                depvar=dependentvar(ocMP(ii));
                OCMATINDIF.endvalue{ii}=depvar(OCMATINDIF.fixendcoordinate{ii},end);
            else
                OCMATINDIF.endvalue{ii}=[];
            end
        case 'l'
    end
    if freeendtime(ii)>0
        depvar=dependentvar(ocMP(ii));
            OCMATINDIF.distance{ii}=norm(OCMATINDIF.saddlepoint{ii}(OCMATINDIF.statecostatecoord)-depvar{end}(OCMATINDIF.statecostatecoord,end));
    else
        OCMATINDIF.distance{ii}=[];
    end
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(time(ocMP(ii).ocgtrajectory));
    OCMATINDIF.arccoord{ii}=[1:arcn]+arcoffset;
    arcoffset=arcoffset+arcn;
end
OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution
sol.parameters=[switchtimes{:}];
counter=0;
paroffset=0;
for ii=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    if ii==1
        if ~isempty(switchtimes{ii})
            OCMATINDIF.switchtimecoord{ii}=(1:length(switchtimes{ii}));
        else
            OCMATINDIF.switchtimecoord{ii}=[];
        end
    else
        OCMATINDIF.switchtimecoord{ii}=paroffset+(1:length(switchtimes{ii}));
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
depvar=dependentvar(ocMP(1));

if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATINDIF.exogenousnumberofstates;
else
    OCMATINDIF.exogenousdynamicscoordinate=[];
end

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.linearization=J;
OCMATINDIF.limitset=limSet;
OCMATINDIF.statecoordinate=scoord;
OCMATINDIF.startvalue=depvar{1}(scoord,1);
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;
OCMATINDIF.objectivevaluecoord=[];


OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=size(OCMATINDIF.asymptoticmatrix{1},2)*indifforder;
OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol=odestruct(ocMultiPath(1));
fn={'x','y','arcarg','arcinterval'};
maxodenum=size(sol.y,1);
for ii=2:nummult
    soltmp=odestruct(ocMultiPath(ii));
    maxodenum=max([maxodenum size(soltmp.y,1)]);
    y=sol.y;
    sol.y=zeros(maxodenum,size(y,2));
    sol.y(1:size(y,1),:)=y;
    y=soltmp.y;
    soltmp.y=zeros(maxodenum,size(y,2));
    soltmp.y(1:size(y,1),:)=y;
    soltmp.x=soltmp.x+sol.x(end);
    soltmp.arcposition=soltmp.arcposition+length(sol.x);
    for jj=1:length(fn)
        sol.(fn{jj})=[sol.(fn{jj}) soltmp.(fn{jj})];
    end
end
