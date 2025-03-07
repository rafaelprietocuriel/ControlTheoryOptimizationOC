function sol=initocmat_AEMF_IS(ocObj,ocMP,targetcoordinate,targetvalue,freesaddlepointcoord,opt,varargin)
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
if nargin==5
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

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=[];
arcoffset=0;
freeparcounter=0;
for ii=1:indifforder
    arcn=arcnum(ocMP(ii));
    solverInfoStruct=solverinfo(ocMP(ii));
    arcarg=arcargument(ocMP(ii));
    if OCMATINDIF.stateconstraint
        sccounter=0;
        arcarg4sc=arcargwithactivestateconstraint(ocObj);
        for jj=2:length(arcarg)
            if ismember(arcarg(jj),arcarg4sc)
                sccounter=sccounter+1;
                OCMATINDIF.jumpcostateindex{ii}(sccounter)=jj;
            end
        end
        if ~sccounter
            OCMATINDIF.jumpcostateindex{ii}=[];
            jump=[];
        end
        if isfield(solverInfoStruct,'jumpcostatecoord')
            OCMATINDIF.jumpcostatecoord{ii}=solverInfoStruct.jumpcostatecoord;
        else
            OCMATINDIF.jumpcostatecoord{ii}=freeparcounter+(1:sccounter);
            freeparcounter=freeparcounter+sccounter;
            sol.parameters=[sol.parameters jump];
        end
    end
    if 0%isfield(solverInfoStruct,'switchtimecoord')
        OCMATINDIF.switchtimecoord{ii}=solverInfoStruct.switchtimecoord;
    else
        arcintv=arcinterval(ocMP(ii));
        OCMATINDIF.switchtimecoord{ii}=freeparcounter+(1:length(arcintv)-2);
        freeparcounter=freeparcounter+length(arcintv)-2;
        sol.parameters=[sol.parameters arcintv(2:end-1)];
    end
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    limSet{ii}=limitset(ocMP(ii));
    J{ii}=linearization(limSet{ii});
    if isempty(inftimetransformation(ocMP(ii)))
        OCMATINDIF.inftimetransformation(ii)=0;
    else
        OCMATINDIF.inftimetransformation(ii)=inftimetransformation(ocMP(ii));
    end
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
    OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(J{ii},pathtpe{ii},'c',ZeroDeviationTolerance);
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcn;
    OCMATINDIF.freesaddlepointcoord{ii}=freesaddlepointcoord{ii};
    OCMATINDIF.equilibriummanifoldcoord{ii}=freeparcounter+(1:length(freesaddlepointcoord{ii}));
    freeparcounter=freeparcounter+length(freesaddlepointcoord{ii});
    if ~isempty(freesaddlepointcoord{ii}) && ~isfield(solverInfoStruct,'freesaddlepointcoord')
        sol.parameters=[sol.parameters OCMATINDIF.saddlepoint{ii}(freesaddlepointcoord{ii}).'];
    end
end

if ~isempty(targetcoordinate)
    sol.parameters=[sol.parameters 0];
    OCMATCONT.continuation=1;
else
    OCMATCONT.continuation=0;
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
if ~isempty(sol.parameters)
    numcontpar=length(continuationparameter(ocMultiPath(1)));
    sol.parameters(end-numcontpar+1:end)=[];
end
x0=initialtime(ocMultiPath);
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
    freepar=parameters(ocMultiPath(ii));
    if ~isempty(freepar)
        numcontpar=length(continuationparameter(ocMultiPath(ii)));
        freepar(end-numcontpar+1:end)=[];
        sol.parameters=[sol.parameters freepar];
    end
end
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
