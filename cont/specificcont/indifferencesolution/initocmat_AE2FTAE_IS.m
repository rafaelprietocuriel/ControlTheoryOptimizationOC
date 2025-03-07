function sol=initocmat_AE2FTAE_IS(ocObj,ocMP,targetcoordinate,targetvalue,opt,varargin)
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
freeendtime=[];

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
targetvalue=targetvalue(:);
stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));

if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if isempty(stopcriterion)
    stopcriterion=0;
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
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(1).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(1).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(1).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(1).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(1).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(1).odedim+(1:domaindata(ii).aedim);
end

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
feidx=zeros(1,indifforder); % index specifying 'finite time equilibrium'
state_costatecoordinate=[statecoord(ocObj);costatecoord(ocObj)];
J=cell(1,indifforder);
pathtpe=cell(1,indifforder);
limSet=cell(1,indifforder);
parameters=[];
for ii=1:indifforder
    depvar=dependentvar(ocMP(ii));
    arcn=arcnum(ocMP(ii));
    arcint=arcinterval(ocMP(ii));
    OCMATINDIF.switchtimecoord{ii}=length(parameters)+(1:arcn-1);
    parameters=[parameters arcint(2:end-1)];
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    if isocasymptotic(ocMP(ii))
        limSet{ii}=limitset(ocMP(ii));
        pathtpe{ii}=pathtype(ocMP(ii));
        if isequilibrium(limSet{ii})
            OCMATINDIF.limitsettype{ii}='e';
            OCMATINDIF.truncationtimecoord{ii}=[];
        elseif isperiodic(limSet{ii})
            OCMATINDIF.limitsettype{ii}='l';
            OCMATINDIF.truncationtimecoord{ii}=solverInfoStruct.truncationtimecoord;
        else
            ocmaterror('Limit set type unknown.')
        end
        J{ii}=linearization(limSet{ii});
        switch OCMATINDIF.limitsettype{ii}
            case 'e'
                OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance);
                OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
            case 'l'
                OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(J{ii},pathtpe{ii},'d',ZeroDeviationTolerance);
                dependentvar(limSet{ii});
                OCMATINDIF.saddlepoint{ii}=ans(:,1);
        end
    else
        feidx(ii)=1;
        OCMATINDIF.asymptoticmatrix{ii}=[];
        OCMATINDIF.saddlepoint{ii}=depvar(state_costatecoordinate,end);
    end
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcoffset+arcn;
end
OCMATINDIF.linearization=J;
OCMATINDIF.ftequilibriumindex=feidx;
OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution

sol.parameters=parameters;

righttimeindex=cumsum(OCMATINDIF.numarc+1);
OCMATINDIF.freeendtime=zeros(1,indifforder);
counter=0;
for ii=1:indifforder
    if feidx(ii)
        OCMATINDIF.truncationtimecoord{ii}=length(sol.parameters)+1;
        sol.parameters=[sol.parameters sol.arcinterval(righttimeindex(ii))];
        OCMATINDIF.freeendtime(ii)=1;
        OCMATINDIF.truncationtime(ii)=NaN;
    elseif ~isempty(freeendtime)
        OCMATINDIF.distance{ii}=norm(OCMATINDIF.saddlepoint{ii}-depvar(:,end));
        OCMATINDIF.truncationtimecoord{ii}=length(sol.parameters)+1;
        sol.parameters=[sol.parameters sol.arcinterval(righttimeindex(ii))];
        OCMATINDIF.freeendtime(ii)=1;
    else
        OCMATINDIF.truncationtime=sol.arcinterval(righttimeindex(ii));
    end
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
end
OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;

OCMATINDIF.statecoordinate=1:statenum(ocObj);
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.startvalue=depvar(targetcoordinate,1);
OCMATINDIF.continuationvector=targetvalue-OCMATINDIF.startvalue;
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;
if ~isempty(targetcoordinate)
    sol.parameters=[sol.parameters 0];
    OCMATCONT.continuation=1;
else
    OCMATCONT.continuation=0;
end

function sol=generatesolstruct(ocMultiPath)

nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
sol.parameters=[];
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
end
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
