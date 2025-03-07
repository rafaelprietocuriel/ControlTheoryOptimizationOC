function sol=initocmat_DAE_IS(docObj,docMP,opt,varargin)
%
% INITOCMAT_DAE_IS initialization for the continuation of an indifference
% threshold
%
% SOL=INITOCMAT_DAE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE) the
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
% SOL=INITOCMAT_DAE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

clear global OCMATCONT OCMATINDIF

global OCMATCONT OCMATINDIF
sol=[];

% input argument docMP is either a cell of docasymptotics or a multi path object 
docMP=docmultipath(docMP);
indifforder=multiplicity(docMP);
targetvalue=[];
targetcoordinate=[];
freevector=[];
targetvectorcoordinate=[];
for ii=1:indifforder
    pathtpe{ii}=pathtype(docMP(ii));
end
if isempty(docObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(docMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
targetvectorcoordinateidx=find(strcmpi(varargin,'targetvectorcoordinate'));
freevectoridx=find(strcmpi(varargin,'freevector'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
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

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     docMP(ii).solver='';
%     docMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(docObj);
OCMATCONT.modelfunc=modelspecificfunc(docObj,'4IndifferenceSolutionContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATINDIF.canonicalsystemmap=funch{1};
OCMATINDIF.canonicalsystemmapjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemmapparameterjacobian=funch{2}{2};
OCMATINDIF.canonicalsystemmaphessian=funch{3}{1};
OCMATINDIF.canonicalsystemmapparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.bcindifference=funch{5}{5};

% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcjacobianasymptotic=funch{6}{2};
OCMATINDIF.bcjacobiantransversality=funch{6}{3};

OCMATINDIF.objectivefunction=funch{8}{1};
OCMATINDIF.objectivefunctionjacobian=funch{8}{2};

% general function
OCMATINDIF.findarcposition=funch{10};
OCMATINDIF.plotcontinuation=funch{11};
OCMATINDIF.testadmissibility=funch{12};
OCMATINDIF.datapath=funch{20};
OCMATINDIF.saveintermediatefiles=funch{21};

sol=generateodiffestruct(docMP);
% mode and path specific variables
OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.targetvectorcoordinate=targetvectorcoordinate;

limSet=cell(1,indifforder);
J=cell(1,indifforder);
OCMATINDIF.arcarg=[];
for ii=1:indifforder
    limSet{ii}=limitset(docMP(ii));
    J{ii}=explicitjacobian(limSet{ii});
    OCMATINDIF.linearization{ii}=linearization(docMP(ii));
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
    OCMATINDIF.asymptoticmatrix{ii}=asymptoticbc(J{ii},pathtpe{ii},'d',ZeroDeviationTolerance);
    OCMATINDIF.saddlepointcodimension(ii)=size(OCMATINDIF.asymptoticmatrix{ii},2);
    OCMATINDIF.arcarg=[OCMATINDIF.arcarg arcargument(docMP(ii))];
    OCMATINDIF.arcnum(ii)=arcnum(docMP(ii));
end
leftside=[1 cumsum(OCMATINDIF.arcnum(1:end-1))+1];
rightside=[leftside(2:end)-1 sum(OCMATINDIF.arcnum)];
OCMATINDIF.arcargcoord=[leftside;rightside];
OCMATINDIF.arcposition=sol.arcposition;
OCMATINDIF.pathcoord=sol.pathposition;

OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.parametervalue=parametervalue(docObj);
OCMATINDIF.initialtime=sol.x(OCMATINDIF.pathcoord(1,:));
OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.statecoordinate=1:statenum(docObj);
OCMATINDIF.statenum=statenum(docObj);
OCMATINDIF.startvalue=sol.y(OCMATINDIF.statecoordinate,1);
OCMATINDIF.continuationvector=targetvalue-OCMATINDIF.startvalue;
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.codimension=1;
OCMATCONT.multipointbvp=true;
OCMATCONT.sumconstraint=true;
OCMATCONT.numsumconstraint=indifforder-1;
OCMATCONT.multipointorder=indifforder;
OCMATCONT.numarc=OCMATINDIF.arcnum;
OCMATCONT.HE.numinitialcondition=numel(targetcoordinate);

function sol=generateodiffestruct(docMP)

sol.x=[];
sol.x0=[];
sol.y=[];
sol.parameters=[];
sol.arcarg=[];
sol.arcposition=[];
rpathposition=[];
for ii=1:multiplicity(docMP)
    sol.x=[sol.x docMP(ii).x0 docMP(ii).x];
    sol.y=[sol.y docMP(ii).y0 docMP(ii).y];
    % add continuation parameter value
    %sol.parameters=[sol.parameters docMP(ii).solverinfo.parameters(1:end-1)];
    sol.arcarg=[sol.arcarg arcargument(docMP(ii))];
    if ii==1
        arcposabs=0;
    else
        arcposabs=arcpos(2,end);
    end
    arcpos=arcposition(docMP(ii));
    arcpos(1,2:end)=arcpos(1,2:end)-1;
    arcpos=arcposabs+arcpos;
    arcpos(2,end)=arcpos(2,end)+1;
    sol.arcposition=[sol.arcposition arcpos];
    if ii<multiplicity(docMP)
        rpathposition=[rpathposition length(docMP(ii).x)+1];
    end
end
lpathposition=[1 rpathposition+1];
rpathposition=[rpathposition length(sol.x)];
sol.pathposition=[lpathposition;rpathposition];

