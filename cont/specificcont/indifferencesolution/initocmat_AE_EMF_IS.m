function sol=initocmat_AE_EMF_IS(ocObj,ocMP,opt,varargin)
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
targetvalue=[];
stopcriterion=[];
freeendtime=[];
targetcoordinate=[];
freevector=[];
targetvectorcoordinate=[];
asymptoticmatrix=[];
explicitemfcoordinate=[];
independentmfcoordinate=[];

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
    opt=defaultocoptions;
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
targetvectorcoordinateidx=find(strcmpi(varargin,'targetvectorcoordinate'));
freevectoridx=find(strcmpi(varargin,'freevector'));
stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'));
explicitemfcoordinateidx=find(strcmpi(varargin,'explicitcoordinate'),1);
indepentmfcoordinateidx=find(strcmpi(varargin,'independentmfcoordinate'),1);

if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
    stopcriterion=0; % stop if limitpoint occurs
end
if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end
if ~isempty(targetvectorcoordinateidx)
    targetvectorcoordinate=varargin{targetvectorcoordinateidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if isempty(stopcriterion)
    stopcriterion=1; % stop if limitpoint occurs
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder); % stop if limitpoint occurs
end
if isempty(asymptoticmatrix)
    asymptoticmatrix=cell(1,indifforder); 
end
if ~isempty(explicitemfcoordinateidx)
    explicitemfcoordinate=varargin{explicitemfcoordinateidx+1};
end
if ~isempty(indepentmfcoordinateidx)
    independentmfcoordinate=varargin{indepentmfcoordinateidx+1};
end

limSet=cell(1,indifforder);
switchtimes=limSet;
limSetDim=zeros(1,indifforder);
if isempty(explicitemfcoordinate)
    explicitemfcoordinate=cell(1,indifforder);
else
    if ~iscell(explicitemfcoordinate) && length(explicitemfcoordinate)==indifforder
        explicitemfcoordinate=num2cell(explicitemfcoordinate);
    elseif length(explicitemfcoordinate)==1
        explicitemfcoordinate=repmat(explicitemfcoordinate,1,indifforder);
    end
end
for ii=1:indifforder
    limSet{ii}=limitset(ocMP(ii));
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
    limSetDim(ii)=length(OCMATINDIF.saddlepoint{ii});
end
if isempty(independentmfcoordinate)
    for ii=1:indifforder
        independentmfcoordinate{ii}=setdiff(1:limSetDim(ii),explicitemfcoordinate{ii});
    end
else
    if ~iscell(independentmfcoordinate) || ~length(independentmfcoordinate)==indifforder
        if ~iscell(independentmfcoordinate) && length(independentmfcoordinate)==indifforder
            independentmfcoordinate=num2cell(independentmfcoordinate);
        elseif length(independentmfcoordinate)==1
            independentmfcoordinate=repmat(independentmfcoordinate,1,indifforder);
        else
            ocmatmsg('Wrong syntax for ''independentmfcoordinate''.')
            return
        end
    end
end
OCMATINDIF.stopcriterion=stopcriterion;
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end

ZeroDeviationTolerance=zeros(1,indifforder);
for ii=1:indifforder
    if iscell(opt)
        ZeroDeviationTolerance(ii)=getocoptions(opt{ii},'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
    else
        ZeroDeviationTolerance(ii)=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
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

OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.targetvectorcoordinate=targetvectorcoordinate;

scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=[];
J=cell(1,indifforder);
arcoffset=0;
parameters=[];
parameterindex=0;
for ii=1:indifforder
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);

    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATINDIF.limitsettype{ii}='e';
    J{ii}=linearization(limSet{ii});
    if isempty(asymptoticmatrix{ii})
        [OCMATINDIF.asymptoticmatrix{ii},numstable,numunstable,numcenter,infoStruct]=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance(ii));


        OCMATINDIF.dependentemfcoordinate{ii}=setdiff(1:limSetDim(ii),[explicitemfcoordinate{ii} independentmfcoordinate{ii}]);
        OCMATINDIF.emfcoord{ii}=setdiff(1:limSetDim(ii),explicitemfcoordinate{ii});
        %OCMATINDIF.emfindex{ii}=parameterindex+(1:length(OCMATINDIF.emfcoord{ii}));
        parameterindex=length(parameters);
        parameters=[parameters OCMATINDIF.saddlepoint{ii}(OCMATINDIF.emfcoord{ii}).'];
        OCMATINDIF.emfindex{ii}=parameterindex+1:length(parameters);
        OCMATINDIF.subspacedim(ii)=numstable;
        OCMATINDIF.orthspacedim(ii)=limSetDim(ii)-OCMATINDIF.subspacedim(ii);
        Y=zeros(OCMATINDIF.orthspacedim(ii),OCMATINDIF.subspacedim(ii));
        OCMATINDIF.Y{ii}=Y;
        %OCMATINDIF.Yindex{ii}=parameterindex+reshape(limSetDim(ii)-length(explicitemfcoordinate{ii})+(1:OCMATINDIF.orthspacedim(ii)*OCMATINDIF.subspacedim(ii)),OCMATINDIF.orthspacedim(ii),OCMATINDIF.subspacedim(ii));
        parameterindex=length(parameters);
        parameters=[parameters Y(:).'];
        OCMATINDIF.Yindex{ii}=parameterindex+1:length(parameters);
        %         Q=computeBase(JacobianMatrix,~stableflag,OCMATINDIF.subspacedim);
        OCMATINDIF.Id{ii}=eye(OCMATINDIF.orthspacedim(ii));
        OCMATINDIF.numY(ii)=numel(Y);
        OCMATINDIF.Q0{ii}=infoStruct.Q;
    else
        OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix{ii};
    end
    if freeendtime(ii)>0
        depvar=dependentvar(ocMP(ii));
        OCMATINDIF.distance{ii}=norm(OCMATINDIF.saddlepoint{ii}-depvar(:,end));
    else
        OCMATINDIF.distance{ii}=[];
    end
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcoffset+arcn;
end
OCMATINDIF.explicitemfcoordinate=explicitemfcoordinate; % a user defined function has to provide the explicit values manifold equilibrium
OCMATINDIF.independentmfcoordinate=independentmfcoordinate; % coordinates used as independent variables for the equilibrium manifold

OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution
counter=0;
for ii=1:indifforder
    parameterindex=length(parameters);
    parameters=[parameters switchtimes{ii}];
    OCMATINDIF.switchtimecoord{ii}=parameterindex+1:length(parameters);
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
end
righttimeindex=cumsum(OCMATINDIF.numarc+1);
OCMATINDIF.truncationtimecoord=[];
for ii=1:indifforder
    if freeendtime(ii)
        parameterindex=length(parameters);
        parameters=[parameters sol.arcinterval(righttimeindex(ii))];
        OCMATINDIF.truncationtimecoord{ii}=parameterindex+1:length(parameters);
        OCMATINDIF.truncationtime{ii}=[];
    else
        OCMATINDIF.truncationtime{ii}=sol.arcinterval(righttimeindex(ii));
    end
end
depvar=dependentvar(ocMP(1));

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];

parameterindex=length(parameters);
parameters=[parameters zeros(1,size(freevector,2))];
OCMATINDIF.freevectorcoord=parameterindex+1:length(parameters);
sol.parameters=parameters;
OCMATINDIF.freevector=freevector;
OCMATINDIF.linearization=J;
OCMATINDIF.limitset=limSet;
OCMATINDIF.statecoordinate=scoord;
OCMATINDIF.startvalue=depvar(scoord,1);
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;

OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=size(OCMATINDIF.asymptoticmatrix{1},2)*indifforder;

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
