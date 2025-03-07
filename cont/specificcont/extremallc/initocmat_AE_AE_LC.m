function sol=initocmat_AE_AE_LC(ocObj,ocAsym,continuationtype,varargin)
%
% initocmat_AE_PER initialization for asymptotic extremal calculation
%
% SOL=initocmat_AE_PER(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
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
fixinitstate=[];
fixtimeindex=[];
continuationcoordinate=[];
continuationparameter=[];
continuationtime=[];
repeatperiod=[];
targetvalue=[];
pathtpe='';
freeendtime=[];
opt=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc equilibrium is empty.\n')
    return
end
asymptoticmatrixidx=find(strcmpi(varargin,'asymptoticmatrix'),1);
repeatperiodidx=find(strcmpi(varargin,'repeatperiod'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate')); % can be empty if type is initialstate, if empty for type parameter is set to the state coordinates
fixtimeindexidx=find(strcmpi(varargin,'fixtimeindex')); % can be empty if type is initialstate, if empty for type parameter is set to the state coordinates
fixdistanceidx=find(strcmpi(varargin,'fixdistance')); % can be empty if type is initialstate, if empty for type parameter is set to the state coordinates
continuationcoordinateidx=find(strcmpi(varargin,'continuationcoordinate')); % has to be specified if type is parameter otherwise empty
continuationparameteridx=find(strcmpi(varargin,'continuationparameter')); % has to be specified if type is initialstate otherwise empty
continuationtimeidx=find(strcmpi(varargin,'continuationtime')); % has to be specified if type is initialstate otherwise empty
targetvalueidx=find(strcmpi(varargin,'targetvalue')); % has to be specified if type is initialstate, can be empty if type is parameter: 
optionidx=find(strcmpi(varargin,'option'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
pathtpeidx=find(strcmpi(varargin,'pathtype'));

if ~isempty(asymptoticmatrixidx)
    asymptoticmatrix=varargin{asymptoticmatrixidx+1};
end
if ~isempty(repeatperiodidx)
    repeatperiod=varargin{repeatperiodidx+1};
end
if ~isempty(pathtpeidx)
    pathtpe=varargin{pathtpeidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end
if ~isempty(fixtimeindexidx)
    fixtimeindex=varargin{fixtimeindexidx+1};
end
if ~isempty(fixdistanceidx)
    fixdistance=varargin{fixdistanceidx+1};
end
if ~isempty(continuationcoordinateidx)
    continuationcoordinate=varargin{continuationcoordinateidx+1};
end
if ~isempty(continuationparameteridx)
    continuationparameter=varargin{continuationparameteridx+1};
end
if ~isempty(continuationtimeidx)
    continuationtime=varargin{continuationtimeidx+1};
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=0;
end
if isempty(pathtpe)
    pathtpe=pathtype(ocAsym);
end
if isempty(pathtpe)
    pathtpe='cs';
end
if isempty(freeendtime)
    freeendtime=0;
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
        if isempty(fixinitstate)
            OCMATAE.fixinitstate=statecoord(ocObj);
        end
    case 'time'
        continuationcoordinate=[];
        if isempty(fixinitstate)
            OCMATAE.fixinitstate=statecoord(ocObj);
        end
    otherwise
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
asymptoticbcmethod=getocoptions(opt,'OCCONTARG','AsymptoticBCMethod');

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

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.bclimitcycle=funch{5}{6};

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
if ~isautonomous(ocObj)
    OCMATAE.canonicalsystemderivativetime=funch{2}{3};
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATAE.objectivefunctionderivativetime=funch{8}{4};
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
OCMATAE.statecostatecoordinate=1:2*statenum(ocObj);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.fixdistance=fixdistance;

if isdynprimitive(ocAsym)
    switch pathtpe
        case {'s','sc','cs'}
            ocLCC{1}=ocAsym;
        case {'u','uc'}
            ocLCC{1}=reverse(ocAsym);
    end
    ocTrj=concatenate(ocLCC{ones(1,repeatperiod)});
    ocAsym=ocasymptotic(ocTrj,ocAsym);
end
switch pathtpe
    case {'s','sc','cs'}
        OCMATAE.reversetime=0;
    case {'u','uc'}
        OCMATAE.reversetime=1;
end
ocLC=limitset(ocAsym);
depvar=dependentvar(ocAsym);
depvarLC=dependentvar(ocLC);

if isempty(asymptoticmatrix)
    JacobianMatrix=jacobian(ocLC);
    [OCMATAE.asymptoticmatrix,OCMATAE.numstable,OCMATAE.numunstable,OCMATAE.numcenter,infoStruct]=asymptoticbc(JacobianMatrix,pathtpe,'d',ZeroDeviationTolerance,asymptoticbcmethod);
end

sol=generateodestruct(ocAsym,continuationtype);

OCMATAE.continuationvector=[];
switch continuationtype
    case 'initialstate'
        arcintv=arcinterval(ocAsym);
        OCMATAE.TRJ.switchtimes=arcintv(2:end-1);
        OCMATAE.TRJ.endtime=arcintv(end);

        OCMATAE.fixtimeindex=fixtimeindex;
        OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATAE.TRJ.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.TRJ.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.TRJ.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.TRJ.freetimecoordinate;
        parametercoordinate=length(sol.parameters);
%         sol.parameters=OCMATAE.TRJ.switchtimes;
%         OCMATAE.TRJ.switchtimecoordinate=parametercoordinate+(1:length(OCMATAE.TRJ.switchtimes));
%         parametercoordinate=length(sol.parameters);
        if OCMATAE.freeendtime
            sol.parameters=[sol.parameters OCMATAE.TRJ.endtime];
            OCMATAE.TRJ.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end

        OCMATAE.TRJ.arcarg=arcargument(ocAsym);
        OCMATAE.TRJ.edge=[OCMATAE.TRJ.arcarg(1:end-1);OCMATAE.TRJ.arcarg(2:end)];

        OCMATAE.TRJ.numarc=arcnum(ocAsym);
        OCMATAE.TRJ.arccoordinate=1:OCMATAE.TRJ.numarc;
        OCMATAE.trajectoryindex=OCMATAE.TRJ.arccoordinate;
        
        %OCMATAE.LC=ocLC;

        OCMATAE.TRJ.continuationcoordinate=continuationcoordinate;
        OCMATAE.startvalue=depvar(continuationcoordinate,1);
        OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point

        OCMATAE.LC.initialpoint=depvarLC(:,1);
        OCMATAE.solutionindex=ones(1,OCMATAE.TRJ.numarc);% 1 ... trajectory, 2 ... limitcycle
        
        % free 

        sol.parameters=[sol.parameters 0];
        OCMATAE.continuationcoordinate=parametercoordinate+1;
        OCMATAE.limitset=ocLC;
        
        OCMATCONT.monodromy=0;

    case 'parameter'
        dxdt=canonicalsystem(ocObj,ocLC,[],1);

        % for the phase condition
        OCMATAE.LC.velocityvector=dxdt(:,1);
        OCMATAE.LC.velocitycoordinate=1:length(dxdt(:,1));
        OCMATAE.LC.velocityvector=OCMATAE.LC.velocityvector/norm(OCMATAE.LC.velocityvector);
        OCMATAE.TRJ.initialpoint=depvar(:,1);
        OCMATAE.LC.initialpoint=depvarLC(:,1);

        arcintv=arcinterval(ocAsym);
        OCMATAE.TRJ.switchtimes=arcintv(2:end-1);
        OCMATAE.TRJ.endtime=arcintv(end);
        arcintv=arcinterval(ocLC);
        OCMATAE.LC.switchtimes=arcintv(2:end-1);
        OCMATAE.LC.period=arcintv(end);

        OCMATAE.fixtimeindex=fixtimeindex;
        OCMATAE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATAE.TRJ.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.TRJ.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.TRJ.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.TRJ.freetimecoordinate;
        parametercoordinate=length(sol.parameters);
%         sol.parameters=OCMATAE.TRJ.switchtimes;
%         OCMATAE.TRJ.switchtimecoordinate=parametercoordinate+(1:length(OCMATAE.TRJ.switchtimes));
%         parametercoordinate=length(sol.parameters);
        if OCMATAE.freeendtime
            sol.parameters=[sol.parameters OCMATAE.TRJ.endtime];
            OCMATAE.TRJ.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end
        sol.parameters=[sol.parameters OCMATAE.LC.switchtimes];
        OCMATAE.LC.switchtimecoordinate=parametercoordinate+(1:length(OCMATAE.LC.switchtimes));
        parametercoordinate=length(sol.parameters);
        sol.parameters=[sol.parameters OCMATAE.LC.period];
        OCMATAE.LC.periodcoordinate=parametercoordinate+1;
        parametercoordinate=length(sol.parameters);


        OCMATAE.TRJ.arcarg=arcargument(ocAsym);
        OCMATAE.TRJ.edge=[OCMATAE.TRJ.arcarg(1:end-1);OCMATAE.TRJ.arcarg(2:end)];

        OCMATAE.LC.arcarg=arcargument(ocLC);
        OCMATAE.LC.edge=[OCMATAE.LC.arcarg(1:end-1);OCMATAE.LC.arcarg(2:end)];

        OCMATAE.TRJ.numarc=arcnum(ocAsym);
        OCMATAE.TRJ.arccoordinate=1:OCMATAE.TRJ.numarc;
        
        OCMATAE.LC.numarc=arcnum(ocLC);
        OCMATAE.LC.arccoordinate=1:OCMATAE.LC.numarc;

        OCMATAE.solutionindex=[ones(1,OCMATAE.TRJ.numarc) ones(1,OCMATAE.LC.numarc)+1];% 1 ... trajectory, 2 ... limitcycle
        OCMATAE.trajectoryindex=OCMATAE.TRJ.arccoordinate;
        OCMATAE.limitcycleindex=OCMATAE.TRJ.arccoordinate(end)+OCMATAE.LC.arccoordinate;
        %OCMATCONT.equilibriumcoordinate=1:canonicalsystemdimension(ocObj);
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

        OCMATAE.continuationparameterindex=parameterindex(ocObj,continuationparameter);
        OCMATAE.modelparameter=modelparameter(ocAsym);
        if ~isempty(targetvalue)
            OCMATAE.startvalue=OCMATAE.modelparameter(OCMATAE.continuationparameterindex);
            OCMATAE.continuationvector=targetvalue(:).'-OCMATAE.startvalue; % continue solution along the line from starting to target value
            sol.parameters=[sol.parameters 0];
        else
            sol.parameters=[sol.parameters OCMATAE.modelparameter(OCMATAE.continuationparameterindex)];
        end
        OCMATAE.continuationcoordinate=parametercoordinate+1;
        OCMATCONT.monodromy=1;
        OCMATCONT.monodromymatrix=JacobianMatrix;
        OCMATCONT.limitcycleposition=sol.arcposition(:,OCMATAE.limitcycleindex);
        OCMATCONT.limitcyclediscretizationlength=OCMATCONT.limitcycleposition(2,end)-OCMATCONT.limitcycleposition(1,1)+1;
    case 'time'
        arcintv=arcinterval(ocAsym);
        OCMATAE.TRJ.switchtimes=arcintv(2:end-1);
        OCMATAE.TRJ.endtime=arcintv(end);

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
        sol.parameters=OCMATAE.TRJ.switchtimes(OCMATAE.freetimeindex);
        OCMATAE.TRJ.freetimecoordinate=parametercoordinate+(1:length(OCMATAE.freetimeindex));
        OCMATAE.TRJ.switchtimecoordinate(OCMATAE.freetimeindex)=OCMATAE.TRJ.freetimecoordinate; % coordinate of the switching times within the free parameters for each arc, 0 ... switching time is fixed
        parametercoordinate=length(sol.parameters);
        if freeendtime
            sol.parameters=[sol.parameters OCMATAE.TRJ.endtime];
            OCMATAE.TRJ.endtimecoordinate=parametercoordinate+1;
            parametercoordinate=length(sol.parameters);
        end

        OCMATAE.TRJ.initialpoint=depvar(:,1);
        OCMATAE.TRJ.arcarg=arcargument(ocAsym);
        OCMATAE.TRJ.edge=[OCMATAE.TRJ.arcarg(1:end-1);OCMATAE.TRJ.arcarg(2:end)];

        OCMATAE.TRJ.numarc=arcnum(ocAsym);
        OCMATAE.TRJ.arccoordinate=1:OCMATAE.TRJ.numarc;
        OCMATAE.trajectoryindex=OCMATAE.TRJ.arccoordinate;

        OCMATAE.LC.initialpoint=depvarLC(:,1);
        OCMATAE.solutionindex=ones(1,OCMATAE.TRJ.numarc);% 1 ... trajectory, 2 ... limitcycle

        % free
        if continuationtime==length(arcintv)-1
            sol.parameters=[sol.parameters OCMATAE.TRJ.endtime];
            OCMATAE.TRJ.endtimecoordinate=parametercoordinate+1;
        else
            sol.parameters=[sol.parameters OCMATAE.TRJ.switchtimes(continuationtime)];
        end
        OCMATAE.continuationtimeindex=continuationtime;
        OCMATAE.continuationcoordinate=parametercoordinate+1;
        OCMATAE.limitset=ocLC;

        OCMATCONT.monodromy=0;

        
end

if OCMATAE.fixdistance
    OCMATAE.distance=norm(OCMATAE.LC.initialpoint-depvar(:,end));
end

% mode and path specific variables
OCMATAE.pathtype=pathtpe;

%OCMATAE.implicitcontrolindex=domaindata(ocPerarcindex).implicitcontrolindex;
pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.codimension=1;
OCBVP.numode=OCMATAE.statecostatecoordinate(end);
% add continuation parameter value

function sol=generateodestruct(ocAsym,continuationtype)
sol=[];
if isempty(ocAsym)
    return
end

switch continuationtype
    case {'time','initialstate'}
        sol.x=independentvar(ocAsym);
        sol.y=dependentvar(ocAsym);
        sol.arcinterval=arcinterval(ocAsym); %
        sol.arcarg=arcargument(ocAsym);
        sol.timehorizon=inf;
        x0=0;
    case 'parameter'
        ocLC=limitset(ocAsym);
        sol.x=independentvar(ocAsym);
        sol.x=[sol.x sol.x(end)+independentvar(ocLC)];
        sol.y=[dependentvar(ocAsym) dependentvar(ocLC)];
        sol.arcinterval=[arcinterval(ocAsym) arcinterval(ocLC)]; %
        sol.arcarg=[arcargument(ocAsym) arcargument(ocLC)];
        sol.timehorizon=[inf period(ocLC)];
        x0=[0 0];
        
    otherwise
end
arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];

