function sol=initocmat_AE_AE_INF(ocObj,ocAsym,continuationtype,varargin)
%
% initocmat_AE_AE_INF initialization for asymptotic extremal calculation
%
% SOL=initocmat_AE_AE_INF(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_AE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
asymptoticmatrix=[];
fixdistance=[];
fixinitstate=[];
fixendcoordinate=[];
fixtimeindex=[];
continuationcoordinate=[];
continuationparameter=[];
continuationtime=[];
targetvalue=[];
pathtpe='';
freeendtime=[];
opt=[];
exogenousfunction=[];
exogenousinitialstates=[];
divergingcoordinate=[];
userbc=[];
objectivevaluecalc=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end

for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        if strcmp(varargin{ii},'option')
            eval('opt=varargin{ii+1};')
        elseif strcmp(varargin{ii},'pathtype')
            eval('pathtpe=varargin{ii+1};')
        else
            eval([varargin{ii} '=varargin{ii+1};'])
        end
    end
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
    pathtpe='s';
end
if isempty(freeendtime)
    freeendtime=0;
end
if isempty(userbc)
    userbc=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end

targetvalue=targetvalue(:);

OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if ~isempty(exogenousfunction)
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
    case 'endstate'
        if isempty(continuationcoordinate)
            ocmatmessage('No continuation coordinates provided.\n')
            return
        end
        if isempty(targetvalue)
            ocmatmessage('No target states provided.\n')
            return
        end
        continuationparameter=[];
        OCMATAE.fixinitstate=fixinitstate;
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
        OCMATAE.fixinitstate=fixinitstate;
    otherwise
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

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
OCMATAE.bcend=funch{5}{5};
if ~isempty(userbc) || userbc
    OCMATAE.userfuncbc=funch{5}{7};
end
OCMATAE.bcinf=funch{5}{8};

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
OCMATAE.statenum=statenum(ocObj);
OCMATAE.statecoordinate=statecoord(ocObj);
OCMATAE.statecostatecoordinate=1:2*statenum(ocObj);
OCMATAE.divergingcoordinate=divergingcoordinate;
OCMATAE.convergingcoordinate=setdiff(OCMATAE.statecostatecoordinate,divergingcoordinate);

OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);
OCMATAE.fixdistance=fixdistance;

OCMATAE.userbc=userbc;

switch pathtpe
    case {'s','sc','cs'}
        OCMATAE.reversetime=0;
    case {'u','uc'}
        OCMATAE.reversetime=1;
end
ocEP=limitset(ocAsym);
depvar=dependentvar(ocAsym);
depvarEP=dependentvar(ocEP);

if isempty(asymptoticmatrix)
    JacobianMatrix=jacobian(ocEP);
    OCMATAE.EP.linearization=JacobianMatrix;
end
Jred=JacobianMatrix;
if ~isempty(Jred)
    Jred(divergingcoordinate,:)=[];
    Jred(:,divergingcoordinate)=[];

    [OCMATAE.asymptoticmatrix,OCMATAE.numstable,OCMATAE.numunstable,OCMATAE.numcenter,infoStruct]=asymptoticbc(Jred,pathtpe,'c',ZeroDeviationTolerance);
else
    OCMATAE.asymptoticmatrix=[];
end
sol=generateodestruct(ocAsym);

numberofodes=canonicalsystemdimension(ocObj);
OCMATAE.objectivevaluecalc=objectivevaluecalc;
if objectivevaluecalc
    o=objectivefunction(ocObj,sol);
    numberofodes=numberofodes+1;
    if size(sol.y,1)<numberofodes
        sol.y(numberofodes,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
    end
    OCMATAE.objectivevaluecoordinate=numberofodes;
else
    if size(sol.y,1)>=numberofodes+1
        sol.y(numberofodes+1,:)=[];
    end    
end
if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATAE.exogenousnumberofstates;
    if size(sol.y,1)<numberofodes+OCMATAE.exogenousnumberofstates
        sol.y(OCMATAE.exogenousdynamicscoordinate,:)=repmat(OCMATAE.exogenousinitialstates,1,size(sol.y,2));
    end

else
    OCMATAE.exogenousdynamicscoordinate=[];
end

OCMATAE.continuationvector=[];
OCMATAE.EP.saddlepoint=depvarEP;
OCMATAE.EP.arcarg=arcargument(ocEP);
OCMATAE.initialpoint=depvar(:,1);
OCMATAE.endpoint=depvar(:,end);
OCMATAE.fixendcoordinate=fixendcoordinate;
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

    case 'endstate' % state or costate
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
        OCMATAE.startvalue=depvar(continuationcoordinate,end);
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
        sol.parameters=[sol.parameters OCMATAE.EP.saddlepoint];
        OCMATAE.EP.saddlepointcoordinate=parametercoordinate+(1:length(OCMATAE.EP.saddlepoint));
        parametercoordinate=length(sol.parameters);


        OCMATAE.arcarg=arcargument(ocAsym);
        OCMATAE.edge=[OCMATAE.arcarg(1:end-1);OCMATAE.arcarg(2:end)];


        OCMATAE.numarc=arcnum(ocAsym);
        OCMATAE.arccoordinate=1:OCMATAE.numarc;

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
        OCMATAE.fixendcoordinate=fixendcoordinate;
        if ~isempty(fixendcoordinate)
            OCMATAE.fixendvalue=depvar(OCMATAE.fixendcoordinate,end);
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
    OCMATAE.distance=norm(OCMATAE.EP.saddlepoint(OCMATAE.convergingcoordinate)-depvar(OCMATAE.convergingcoordinate,end));
end

% mode and path specific variables
OCMATAE.pathtype=pathtpe;

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);


OCMATCONT.codimension=1;
OCBVP.numode=OCMATAE.statecostatecoordinate(end);

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

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           