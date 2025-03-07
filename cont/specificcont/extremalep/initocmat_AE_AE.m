function sol=initocmat_AE_AE(ocObj,ocAsym,continuationtype,varargin)
%
% initocmat_AE_AE initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_AE(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
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
%
% SHORT EXPLANATION OF THE MATHEMATICAL BACKGROUND
%
% The underlying mathematical problem formulation is to find a trajectory
% (time path) x(t)=(x_1(t),...,x_N(t)) with initial condition x_j(0)=x_j^0,
% j=i_1,...,i_k and lim_{t\to\infty}x(t)=hat-x (convergence to the
% equilibrium).
%
% The problem is solved using a continuation algorithm, where the
% continuation is done for the initial condition
%   x_j(0)=x_j^0, j=i_1,...,i_k
% and the continuation parameter 'mu' is defined as
%   x_j(0)=x_j^0*mu+(1-mu)*hat-x_j, j=i_1,...,i_k.
% Thus, for mu=0 we have
%   x_j(0)=hat-x_j, j=i_1,...,i_k
% and for mu=1
%   x_j(0)=x_j^0
% The end condition, convergence to the equilibrium, is reformulated in a
% way that allows a numerical treatment. The default way is the truncation
% of the infinite time to a finite time 'T' and the condtion that the end
% point x(T) ends on the linearized stable manifold (stable eigenspace).
%
% This means that at the start of the continuation the equilibrium path
% (constant solution at the equilibrium) trivially satisfies the initial
% and end condition. Therefore, with the provision of the equilibrium the
% initial solution is given as well. In that sense the OCEP argument
% performs two tasks. The searched for solution converges to OCEP and OCEP
% is the initial solution (mu=0) of the continuation process.
%
% The denomination as TARGETVALUE maybe misleading but from the
% continuation point of view it denotes the target. From the problem
% perspective it denotes the initial state values of the searched for
% solution.
%
% SOL=INITOCMAT_AE_EP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE,OPT) with the
% option structure OPT the threshold 'ZeroDeviationTolerance' and initial
% number of discretization points 'TrivialArcMeshNum' for the equilibrium
% solution can be changed.
%   ZeroDeviationTolerance ... provides the tolerance to classify an
%                              eigenvalue numerically as zero.
%   TrivialArcMeshNum      ... provides the number of points for the
%                              constant solution at the equilbrium.
%
% SOL=INITOCMAT_AE_EP(...,'TruncationTime',T) the truncation of the
% infinite time horizon to the finite time T
%
% SOL=INITOCMAT_AE_EP(...,'PathType',p)
%   p='s' (default) stable saddle-path calculation
%   p='u' unstable saddle-path calculation

clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];

asymptoticmatrix=[];
fixdistance=[];
fixdistancecoordinate=[];
fixinitstate=[];
fixtimeindex=[];
continuationcoordinate=[];
continuationparameter=[];
targetcontrolcoordinate=[];
targetcontrolvalue=[];
continuationtime=[];
truncationtime=[];
targetvalue=[];
pathtpe='';
freeendtime=[];
option=[];
objectivevaluecalc=[];
exogenousfunction=[];
exogenousinitialstates=[];
userbc=[];
userbcvalue=[];
hitfunction=[];
targettype='';
freeparameter=[];
simple=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
OCMATAE.implicit=false;
if isimplicit(ocObj)
    OCMATAE.implicit=true;
    sol=initocmat_AE_AE4implicit(ocObj,ocAsym,continuationtype,varargin{:});
    return
else
    OCMATAE.EP.implicitcontrolnum=0;
end
% call all provided specifications
for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        if strcmp(varargin{ii},'pathtype')
            pathtpe=varargin{ii+1};
        else
            eval([varargin{ii} '=varargin{ii+1};'])
        end
    end
end

if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=0;
end
if isempty(pathtpe)
    pathtpe=pathtype(ocAsym);
end
if isempty(pathtpe)
    pathtpe='s';
end
if isempty(freeendtime)
    freeendtime=0;
end
if isempty(userbc)
    userbc=0;
end
if isempty(hitfunction)
    hitfunction=0;
end
OCMATAE.hitfunction=hitfunction;
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(simple)
    simple=0;
end
targetvalue=targetvalue(:);

OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if OCMATAE.exogenousfunction
    OCMATAE.exogenousinitialstates=exogenousinitialstates;
    OCMATAE.exogenousnumberofstates=length(exogenousinitialstates);
end

OCMATCONT.continuationtype=continuationtype;
switch continuationtype
    case 'initialstate'
        if isempty(continuationcoordinate)
            ocmatmessage('No continuation coordinates provided.\n')
            return
        end
        if isempty(targetvalue)
            ocmatmessage('No target states provided.\n')
            return
        end
        continuationparameter=[];
    case 'parameter'
        if isempty(continuationparameter)
            ocmatmessage('No continuation parameter provided.\n')
            return
        end
        continuationcoordinate=[];
        if isempty(fixinitstate) && ~userbc
            OCMATAE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATAE.initialstatecoordinate=fixinitstate;
        end
    case 'time'
        continuationcoordinate=[];
        if isempty(fixinitstate) && ~userbc
            OCMATAE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATAE.initialstatecoordinate=fixinitstate;
        end
        if isempty(targettype)
            targettype='T';
        end
        OCMATAE.targettype=targettype;
    otherwise
end
OCMATCONT.ZeroDeviationTolerance=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
TrivialArcMeshNum=getocoptions(option,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)
OCMATCONT.AsymptoticBCMethod=getocoptions(option,'OCCONTARG','AsymptoticBCMethod');
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};

if hitfunction
    OCMATAE.targetfunction=funch{4};
end

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.equilibrium=funch{5}{4};
if ~isempty(userbc) && userbc
    OCMATAE.userbcfunc=funch{5}{7};
end

if ~isempty(targetcontrolcoordinate)
    OCMATAE.bccontrol=funch{5}{8};
end
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
if objectivevaluecalc
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATAE.objectivefunctionderivativetime=funch{8}{4};
end
if ~isautonomous(ocObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
end
if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamics=funch{4}{1};
    OCMATAE.exogenousjacobian=funch{4}{2};
    OCMATAE.exogenousparameterjacobian=funch{4}{3};
end
% function for workspace initialization
try
    OCMATAE.probleminit=funch{10};
catch
    OCMATAE.probleminit=[];
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

OCMATAE.freeendtime=freeendtime;
targetvalue=targetvalue(:);
OCMATAE.targetvalue=targetvalue;

% model information
OCMATAE.statecoordinate=statecoord(ocObj);
coord=costatecoordinate(ocObj);
OCMATAE.statecostatecoordinate=1:coord(end);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.fixdistance=fixdistance;
if fixdistance && isempty(fixdistancecoordinate)
    fixdistancecoordinate=OCMATAE.statecostatecoordinate;
end
OCMATAE.fixdistancecoordinate=fixdistancecoordinate;
OCMATAE.simple=simple;

OCMATAE.userbc=userbc;

switch pathtpe
    case {'s','sc','cs'}
        OCMATAE.reversetime=0;
    case {'u','uc','stu'}
        OCMATAE.reversetime=1;
        truncationtime=-truncationtime;
end

if isdynprimitive(ocAsym)
    if isequilibrium(ocAsym)
        if isempty(truncationtime)
            return
        end
        ocTrj=equilibrium2trajectory(ocAsym,TrivialArcMeshNum,truncationtime);
        ocAsym=ocasymptotic(ocTrj,ocAsym);
    else
        return
    end
end
ocEP=limitset(ocAsym);
depvar=dependentvar(ocAsym);
depvarEP=dependentvar(ocEP);

if isempty(asymptoticmatrix)
    JacobianMatrix=jacobian(ocEP);
    [OCMATAE.asymptoticmatrix,OCMATAE.numstable,OCMATAE.numunstable,OCMATAE.numcenter,infoStruct]=asymptoticbc(JacobianMatrix,pathtpe,'c',OCMATCONT.ZeroDeviationTolerance,OCMATCONT.AsymptoticBCMethod);
    OCMATAE.EP.linearization=JacobianMatrix;
end

sol=generateodestruct(ocAsym);

numberofodes=canonicalsystemdimension(ocObj);
if objectivevaluecalc && ~isfield(ocAsym.solverinfo,'objectivevaluecoordinate')
    o=objectivefunction(ocObj,ocAsym,1);
    OCMATAE.objectivevaluecoordinate=numberofodes;
    if isa(ocObj,'differentialgame')
        for jj=1:2
            sol.y(end+1,:)=[0 cumsum((o(jj,1:end-1)+o(jj,2:end))/2.*diff(time(ocObj,sol,1)))];
            numberofodes=numberofodes+1;
        end
    else
        o=objectivefunction(ocObj,sol,1);
        sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
        numberofodes=numberofodes+1;
    end
    OCMATAE.objectivevaluecoordinate=OCMATAE.objectivevaluecoordinate+1:numberofodes;
    OCMATAE.objectivevaluecoordinate=OCMATAE.objectivevaluecoordinate(:);
elseif isfield(ocAsym.solverinfo,'objectivevaluecoordinate')
    OCMATAE.objectivevaluecoordinate=ocAsym.solverinfo.objectivevaluecoordinate;
end

if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATAE.exogenousnumberofstates;
else
    OCMATAE.exogenousdynamicscoordinate=[];
end
parametercoordinate=length(sol.parameters);
if ~isempty(userbcvalue)
    sol.parameters=userbcvalue;
    OCMATAE.userbccoordinaten=parametercoordinate+(1:length(OCMATAE.userbcvalue));
end
OCMATAE.freeparameter=freeparameter;
OCMATAE.fixinitstate=fixinitstate;
OCMATAE.continuationvector=[];
OCMATAE.EP.saddlepoint=depvarEP;
OCMATAE.EP.arcarg=arcargument(ocEP);
OCMATAE.initialpoint=depvar(:,1);

switch continuationtype
    case 'initialstate'
        arcintv=arcinterval(ocAsym);
        OCMATAE.switchtimes=arcintv(2:end-1);
        OCMATAE.endtime=arcintv(end);

        OCMATAE.fixtimeindex=fixtimeindex;
        OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATAE.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.freetimecoordinate;
        parametercoordinate=length(sol.parameters);
        if OCMATAE.freeendtime
            sol.parameters=[sol.parameters OCMATAE.endtime];
            OCMATAE.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end

        OCMATAE.arcarg=arcargument(ocAsym);
        OCMATAE.edge=[OCMATAE.arcarg(1:end-1);OCMATAE.arcarg(2:end)];

        OCMATAE.numarc=arcnum(ocAsym);
        OCMATAE.arccoordinate=1:OCMATAE.numarc;
        OCMATAE.trajectoryindex=OCMATAE.arccoordinate;

        OCMATAE.depvarcoordinate=continuationcoordinate;
        OCMATAE.startvalue=depvar(continuationcoordinate,1);
        OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point

        OCMATAE.EP.saddlepoint=depvarEP;

        % free

        sol.parameters=[sol.parameters 0];
        OCMATAE.continuationcoordinate=parametercoordinate+1;

    case 'parameter'

        arcintv=arcinterval(ocAsym);
        OCMATAE.switchtimes=arcintv(2:end-1);
        OCMATAE.endtime=arcintv(end);

        OCMATAE.fixtimeindex=fixtimeindex;
        OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATAE.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.freetimecoordinate;
        parametercoordinate=length(sol.parameters);
        if OCMATAE.freeendtime
            sol.parameters=[sol.parameters OCMATAE.endtime];
            OCMATAE.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end
        sol.parameters=[sol.parameters OCMATAE.EP.saddlepoint.'];
        OCMATAE.EP.equilibriumcoordinate=parametercoordinate+(1:length(OCMATAE.EP.saddlepoint));
        parametercoordinate=length(sol.parameters);


        OCMATAE.arcarg=arcargument(ocAsym);
        OCMATAE.edge=[OCMATAE.arcarg(1:end-1);OCMATAE.arcarg(2:end)];


        OCMATAE.numarc=arcnum(ocAsym);
        OCMATAE.arccoordinate=1:OCMATAE.numarc;
        OCMATAE.trajectoryindex=OCMATAE.arccoordinate;

        if ~OCMATAE.simple
            switch pathtpe
                case 's'
                    OCMATAE.subspacedim=OCMATAE.numstable;
                case 'u'
                    OCMATAE.subspacedim=OCMATAE.numunstable;
                case {'sc','cs'}
                    OCMATAE.subspacedim=OCMATAE.numstable+OCMATAE.numcenter;
            end
            OCMATAE.orthspacedim=OCMATAE.statecostatecoordinate(end)-OCMATAE.subspacedim;
            Y=zeros(OCMATAE.orthspacedim,OCMATAE.subspacedim);
            OCMATAE.Y=Y;
            % orthogonal basis for stable eigenspace
            OCMATAE.Q0=infoStruct.Q;
            OCMATAE.Id=eye(OCMATAE.orthspacedim);
            OCMATAE.numY=numel(Y);

            %[R11,R12,R21,R22]=RicattiCoeff(OCMATAE.Q0,JacobianMatrix,OCMATAE.subspacedim);
            sol.parameters=[sol.parameters Y(:).'];
            OCMATAE.Ycoordinate=reshape(parametercoordinate+(1:OCMATAE.orthspacedim*OCMATAE.subspacedim),OCMATAE.orthspacedim,OCMATAE.subspacedim);
            parametercoordinate=length(sol.parameters);
        end
        OCMATAE.continuationparameterindex=parameterindex(ocObj,continuationparameter);
        OCMATAE.freeparameterindex=parameterindex(ocObj,continuationparameter);
        OCMATAE.modelparameter=modelparameter(ocAsym);
        if ~isempty(targetvalue)
            OCMATAE.startvalue=OCMATAE.modelparameter(OCMATAE.continuationparameterindex);
            OCMATAE.continuationvector=targetvalue(:).'-OCMATAE.startvalue; % continue solution along the line from starting to target value
            sol.parameters=[sol.parameters 0];
        else
            sol.parameters=[sol.parameters OCMATAE.modelparameter(OCMATAE.freeparameterindex)];
        end
        OCMATAE.freeparametercoordinate=parametercoordinate+(1:length(OCMATAE.freeparameterindex));
        OCMATAE.continuationcoordinate=length(sol.parameters);
    case 'time'
        arcintv=arcinterval(ocAsym);
        OCMATAE.switchtimes=arcintv(2:end-1);
        OCMATAE.endtime=arcintv(end);

        OCMATAE.fixtimeindex=fixtimeindex;

        if isempty(continuationtime)
            continuationtime=length(arcintv)-1;
        end
        if continuationtime==length(arcintv)-1
            if freeendtime
                ocmatmsg('End time is continued therefore ''freeendtime'' is set to 0.\n')
                freeendtime=0;
            end
            OCMATAE.freeendtime=1;
            OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        else
            OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,[fixtimeindex continuationtime]);
        end


        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATAE.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.freetimecoordinate; % coordinate of the switching times within the free parameters for each arc, 0 ... switching time is fixed
        parametercoordinate=length(sol.parameters);
        if freeendtime
            sol.parameters=[sol.parameters OCMATAE.endtime];
            OCMATAE.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end


        OCMATAE.arcarg=arcargument(ocAsym);
        OCMATAE.edge=[OCMATAE.arcarg(1:end-1);OCMATAE.arcarg(2:end)];

        OCMATAE.numarc=arcnum(ocAsym);
        OCMATAE.arccoordinate=1:OCMATAE.numarc;
        OCMATAE.trajectoryindex=OCMATAE.arccoordinate;

        % free
        if continuationtime==length(arcintv)-1
            sol.parameters=[sol.parameters OCMATAE.endtime];
            OCMATAE.endtimecoordinate=parametercoordinate+1;
        else
            sol.parameters=[sol.parameters OCMATAE.switchtimes(continuationtime)];
        end
        OCMATAE.continuationtimeindex=continuationtime;
        OCMATAE.continuationcoordinate=parametercoordinate+1;
end

if OCMATAE.fixdistance
    OCMATAE.distance=norm(OCMATAE.EP.saddlepoint(OCMATAE.fixdistancecoordinate)-depvar(OCMATAE.fixdistancecoordinate,end));
end

OCMATAE.targetcontrolcoordinate=targetcontrolcoordinate;
if ~isempty(targetcontrolcoordinate)
    if isempty(targetcontrolvalue) || length(targetcontrolcoordinate)~=length(targetcontrolvalue)
        ocmatmessage('Dimensions of ''targetcontrolcoordinate'' and ''targetcontrolvalue'' must agree.')
        sol=[];
        return
    end
    ctrl=control(ocObj,ocAsym);
    OCMATAE.initialcontrolvalue=ctrl(OCMATAE.targetcontrolcoordinate,1);
    OCMATAE.controlvaluedirection=targetcontrolvalue(:)-OCMATAE.initialcontrolvalue;
end
% mode and path specific variables
OCMATAE.pathtype=pathtpe;

%OCMATAE.implicitcontrolindex=domaindata(ocPerarcindex).implicitcontrolindex;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATAE.objectivevaluecalc=objectivevaluecalc;

OCMATCONT.codimension=1;
OCBVP.numode=OCMATAE.statecostatecoordinate(end);
% add continuation parameter value

function sol=generateodestruct(ocAsym)
sol=[];
if isempty(ocAsym)
    return
end

sol.x=independentvar(ocAsym);
sol.y=dependentvar(ocAsym);
sol.arcinterval=arcinterval(ocAsym); %
sol.arcarg=arcargument(ocAsym);
sol.timehorizon=inf;
x0=0;

arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];

