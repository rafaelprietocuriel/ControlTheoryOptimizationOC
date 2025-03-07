function sol=initocmat_AE_SPEC(ocObj,ocTrj,contcoordinate,targetvalue,opt,varargin)
%
% initocmat_AE_SPEC initialization for asymptotic extremal calculation with
% user specific boundary conditions
%
% SOL=initocmat_AE_SPEC(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
% saddle path calculation is initialized.
% OCOBJ          ... corresponding optimal control model
% OCEP           ... equilibrium point (dynprimitive) (hat-x) with (local)
%                    stable manifold of dimension k
% CONTCOORDINATE ... coordinates i_1,...,i_k of the continuation variable
%                    (usually state coordinate(s) in optimal control
%                    problems)
% TARGETVALUE    ... determines direction of the continuation (x_j^0,
%                    j=i_1,...,i_k)
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and OCMATAE are
% initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.

clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
freeparameter='';
freeparametervalue=[];
objectivevaluecalc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocTrj)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end

freeparameteridx=find(strcmpi(varargin,'freeparameter'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end

if any(strcmpi(freeparameter,'timehorizon'))
    timehorizon=arcinterval(ocTrj);
    timehorizon=timehorizon(end);
    OCMATAE.freetimehorizon=1;
else
    timehorizon=[];
    OCMATAE.freetimehorizon=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
targetvalue=targetvalue(:);

% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4UserSpecificContinuation');

% initialize global variable (OCMATAE) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemderivativetime=funch{2}{3};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcspecific=funch{5}{1};

% function for Jacobian
OCMATAE.bcjacobianspecific=funch{6}{1};

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
OCMATAE.objectivefunction=funch{8}{1};
OCMATAE.objectivefunctionjacobian=funch{8}{2};
OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specific model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

% initialize solution using the repeated entries of the equilibrium
sol=generateodestruct(ocTrj,timehorizon);
OCMATAE.timehorizoncoord=length(sol.parameters);

% add continuation parameter value
sol.parameters=[sol.parameters 0];
if objectivevaluecalc
    o=objectivevalue(ocObj,sol);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol)))];
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% mode and path specific variables
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.arcinterval)-2;
if objectivevaluecalc
    OCMATAE.objectivevaluecoord=size(sol.y,1);
end
OCMATAE.objectivevaluecalc=objectivevaluecalc;

OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=sol.y(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.codimension=1;

function sol=generateodestruct(ocTrj,T)

sol=[];
if isempty(ocTrj)
    return
end
repmat(0:T-1,length(independentvar(ocTrj)),1);
sol.x=ocTrj.x;
sol.y=ocTrj.y;
sol.parameters=arcinterval(ocTrj);
sol.parameters(1)=[];
sol.x0=initialtime(ocTrj);
sol.arcinterval=arcinterval(ocTrj); %
sol.timehorizon=inf;
sol.arcarg=arcargument(ocTrj);
arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.idata.tangent=[];
sol.idata.coeff=[];