function sol=initocmat_AE_PER_IS(ocObj,ocMP,targetcoordinate,targetvalue,opt,varargin)
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
targetvalue=targetvalue(:);

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
if ~isautonomous(ocObj)
    OCMATINDIF.canonicalsystemderivativetime=funch{2}{3};
    OCMATINDIF.objectivefunction=funch{8}{1};
    OCMATINDIF.objectivefunctionjacobian=funch{8}{2};
    OCMATINDIF.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATINDIF.objectivefunctionderivativetime=funch{8}{4};
end
OCMATINDIF.autonomous=isautonomous(ocObj);
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
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
    OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

r=discountrate(ocObj);

sol=generatesolstruct(ocMP);
if ~isempty(targetcoordinate)
    sol.parameters=[sol.parameters 0];
    OCMATCONT.continuation=1;
else
    OCMATCONT.continuation=0;
end
OCMATINDIF.objectivevaluecoord=size(sol.y,1);
% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
for ii=1:indifforder
    arcn=arcnum(ocMP(ii));
    solverInfoStruct=solverinfo(ocMP(ii));
    if 0%isfield(solverInfoStruct,'switchtimecoord')
        OCMATINDIF.switchtimecoord{ii}=solverInfoStruct.switchtimecoord;
    else
        if ii==1
            OCMATINDIF.switchtimecoord{ii}=1:arcn-1;
        else
            lastcoord=OCMATINDIF.switchtimecoord{ii-1};
            if ~isempty(lastcoord)
                lastcoord=lastcoord(end);
            else
                lastcoord=0;
            end
            OCMATINDIF.switchtimecoord{ii}=lastcoord+(1:arcn-1);
        end
    end
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    limSet{ii}=limitset(ocMP(ii));
    tau=period(limSet{ii});
    ov=objectivevalue(ocObj,limSet{ii});
    ov=ov(end);
    arct=arcinterval(ocMP(ii));
    arct=arct(end);
    J{ii}=jacobian(limSet{ii});
    if isempty(inftimetransformation(ocMP(ii)))
        OCMATINDIF.inftimetransformation(ii)=0;
    else
        OCMATINDIF.inftimetransformation(ii)=inftimetransformation(ocMP(ii));
    end
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
    OCMATINDIF.saddlepoint{ii}=OCMATINDIF.saddlepoint{ii}(:,end);
    OCMATINDIF.limitset{ii}=limSet{ii};
    OCMATINDIF.linearization{ii}=J{ii};
    OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(J{ii},pathtpe{ii},'d',ZeroDeviationTolerance);
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    OCMATINDIF.endcoord(ii)=OCMATINDIF.arccoord{ii}(end);
    arcoffset=OCMATINDIF.arccoord{ii}(end);
end
OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution
counter=0;
for order=1:indifforder
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(order);
    OCMATINDIF.solutionindex(counter_start:counter)=order;
end
depvar=dependentvar(ocMP(1));

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
righttimeindex=cumsum(OCMATINDIF.numarc+1);
OCMATINDIF.truncationtime=sol.arcinterval(righttimeindex);
OCMATINDIF.linearization=J;
OCMATINDIF.statecoordinate=1:statenum(ocObj);
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.startvalue=depvar(targetcoordinate,1);
OCMATINDIF.continuationvector=targetvalue-OCMATINDIF.startvalue;
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=1;

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
numcontpar=length(continuationparameter(ocMultiPath(1)));
sol.parameters(end-numcontpar+1:end)=[];
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end-1);
end
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
    freepar=parameters(ocMultiPath(ii));
    if isempty(freepar)
        freepar=actarcinterval(2:end-1);
        numcontpar=0;
    else
        numcontpar=length(continuationparameter(ocMultiPath(ii)));

    end
    freepar(end-numcontpar+1:end)=[];
    sol.parameters=[sol.parameters freepar];
end
sol.x0=initialtime(ocMultiPath(1));
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
