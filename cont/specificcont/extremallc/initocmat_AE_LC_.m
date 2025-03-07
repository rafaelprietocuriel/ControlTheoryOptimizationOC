function sol=initocmat_AE_LC(ocObj,ocAsym,fixcoordinate,opt,varargin)
%
% INITOCMAT_AE_IS initialization for the continuation of an indifference
% threshold
%
% SOL=INITOCMAT_AE_IS(OCOBJ,ocAsym,TARGETCOORDINATE,TARGETVALUE) the
% continuation of an indifference threshold in the state space is
% initialized.
% OCOBJ          ... corresponding optimal control model
% ocAsym           ... two cell array of ocasymptotics, for the different
%                    solution paths or an instance of an ocAsym
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
% SOL=INITOCMAT_AE_IS(OCOBJ,ocAsym,TARGETCOORDINATE,TARGETVALUE,OPT)

clear global OCMATCONT OCMATLC

global OCMATCONT OCMATLC
sol=[];

% input argument ocAsym is either a cell of ocasymptotics or a multi path object
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
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

if nargin==4
    opt=defaultocoptions;
end


OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATLC) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATLC.canonicalsystem=funch{1};
OCMATLC.canonicalsystemjacobian=funch{2}{1};
OCMATLC.canonicalsystemparameterjacobian=funch{2}{2};
OCMATLC.canonicalsystemhessian=funch{3}{1};
OCMATLC.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATLC.bcinitial=funch{5}{1};
OCMATLC.bcasymptotic=funch{5}{2};
OCMATLC.bctransversality=funch{5}{3};
OCMATLC.bcindifference=funch{5}{5};

% function for Jacobian
OCMATLC.bcjacobianinitial=funch{6}{1};
OCMATLC.bcjacobianasymptotic=funch{6}{2};
OCMATLC.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATLC.hybridinfo=funch{7}{1};
OCMATLC.domain=funch{7}{2};
OCMATLC.guard=funch{7}{3};
OCMATLC.reset=funch{7}{4};
OCMATLC.switchtime=funch{7}{5};
OCMATLC.jacobianguard=funch{7}{7};
OCMATLC.jacobianreset=funch{7}{8};
OCMATLC.domaindiscretization=funch{7}{9};
OCMATLC.timesettransformation=funch{7}{10};

% general function
OCMATLC.plotcontinuation=funch{11};
OCMATLC.testadmissibility=funch{12};
OCMATLC.datapath=funch{20};
OCMATLC.saveintermediatefiles=funch{21};

hybridinfo=OCMATLC.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATLC.domain(hybridinfo.arcarg(ii));
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

sol=generatesolstruct(ocAsym);
% mode and path specific variables

sol.parameters=sol.arcinterval(2:end);

numswitchtimes=length(sol.arcinterval(2:end-1));
OCMATLC.switchtimecoord=1:numswitchtimes;
OCMATLC.periodcoord=numswitchtimes+1;

depvar=dependentvar(ocAsym);

OCMATLC.statecoord=statecoord(ocObj);
OCMATLC.costatecoord=costatecoord(ocObj);
OCMATLC.initialcostatedifference=ocAsym.y(OCMATLC.costatecoord,1)-ocAsym.y(OCMATLC.costatecoord,end);
OCMATLC.initialstatedifference=ocAsym.y(OCMATLC.statecoord,1)-ocAsym.y(OCMATLC.statecoord,end);

OCMATLC.parametervalue=parametervalue(ocObj);
OCMATLC.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATLC.numarc(1:end-1)+1)];
%righttimeindex=cumsum(OCMATLC.numarc+1);
OCMATLC.statecoordinate=1:statenum(ocObj);
OCMATLC.fixcoordinate=fixcoordinate;
OCMATLC.fixvalue=depvar(fixcoordinate,1);
pathname=OCMATLC.datapath();
[resultfile,globalvarfile]=OCMATLC.saveintermediatefiles();
OCMATLC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATLC.objectivevaluecalc=0;
OCMATLC.autonomous=1;

sol.parameters=[sol.parameters 1];

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocAsym)

sol.x=independentvar(ocAsym);
sol.x0=sol.x(1);
sol.y=dependentvar(ocAsym);
sol.arcarg=arcargument(ocAsym);
sol.arcinterval=arcinterval(ocAsym);
%sol.parameters=parameters(ocAsym(1));
%numcontpar=length(continuationparameter(ocAsym(1)));
%sol.parameters(end-numcontpar+1:end)=[];
%if isempty(sol.parameters)
%    sol.parameters=sol.arcinterval(2:end-1);
%end
sol.parameters=[];
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];


