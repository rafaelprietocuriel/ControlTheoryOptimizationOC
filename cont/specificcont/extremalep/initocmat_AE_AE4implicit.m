function sol=initocmat_AE_AE4implicit(ocObj,ocgAsym,continuationtype,varargin)
%

global OCMATCONT OCMATAE
sol=[];

asymptoticmatrix=[];
fixdistance=[];
fixdistancecoordinate=[];
fixinitstate=[];
fixtimeindex=[];
continuationcoordinate=[];
continuationparameter=[];
continuationtime=[];
targetcontrolcoordinate=[];
targetcontrolvalue=[];
truncationtime=[];
targetvalue=[];
pathtpe='';
targettype='';
freeendtime=[];
option=[];
objectivevaluecalc=[];
exogenousfunction=[];
hitfunction=[];
userbc=[];
freeparameter=[];
simple=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end
if isempty(ocgAsym)
    ocmatmsg('Initial solution path is empty.\n')
    return
end

% call all provided specifications
for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end

if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=0;
end
if isempty(targettype)
    targettype='T';
end
if isempty(pathtpe)
    pathtpe=pathtype(ocgAsym);
    if isempty(pathtpe)
        pathtpe='s';
    end
end
if isempty(freeendtime)
    freeendtime=0;
end
if isempty(userbc)
    userbc=0;
end
if isempty(exogenousfunction)
    exogenousfunction=0;
end

if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end

if isempty(hitfunction)
    hitfunction=0;
end

if isempty(simple)
    simple=0;
end
% Exogenous functions
OCMATAE.exogenousfunction=exogenousfunction;
OCMATAE.exogenousnumberofstates=0;
if exogenousfunction
    OCMATAE.exogenousinitialstates=exogenousinitialstates;
    OCMATAE.exogenousnumberofstates=length(exogenousinitialstates);
end
OCMATAE.objectivevaluecalc=objectivevaluecalc;
OCMATAE.targettype=targettype;
OCMATAE.fixinitstate=fixinitstate;
OCMATAE.simple=simple;
OCMATAE.targetcontrolcoordinate=targetcontrolcoordinate;

% set options
OCMATCONT.ZeroDeviationTolerance=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
TrivialArcMeshNum=getocoptions(option,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)
OCMATCONT.AsymptoticBCMethod=getocoptions(option,'OCCONTARG','AsymptoticBCMethod');


OCMATAE.pathtype=pathtpe;
switch OCMATAE.pathtype
    case {'s','sc','cs'}
        OCMATAE.reversetime=0;
    case {'u','uc'}
        OCMATAE.reversetime=1;
        if ~isempty(truncationtime)
            truncationtime=-truncationtime;
        end
end

% 
if isgdynprimitive(ocgAsym)
    if isequilibrium(ocgAsym)
        if isempty(truncationtime)
            return
        end
        ocTrj=equilibrium2trajectory(ocgAsym,TrivialArcMeshNum,truncationtime);
        ocgAsym=ocgasymptotic(ocTrj,ocgAsym);
    else
        return
    end
end
% model information
OCMATAE.statecoordinate=statecoord(ocObj);
OCMATAE.statecostatecoordinate=1:2*statenum(ocObj);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.autonomous=isautonomous(ocObj);

% initialize global variable (OCMATAE) for specific continuation process
OCMATAE.freeendtime=freeendtime;
OCMATAE.targetvalue=targetvalue(:);
OCMATAE.fixdistance=fixdistance;
if fixdistance && isempty(fixdistancecoordinate)
    fixdistancecoordinate=OCMATAE.statecostatecoordinate;
end
OCMATAE.fixdistancecoordinate=fixdistancecoordinate;
OCMATAE.userbc=userbc;
OCMATAE.hitfunction=hitfunction;
OCMATCONT.continuationtype=continuationtype;



% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.canonicalsystem=funch{1};
OCMATAE.canonicalsystemjacobian=funch{2}{1};
OCMATAE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATAE.canonicalsystemhessian=funch{3}{1};
OCMATAE.canonicalsystemparameterhessian=funch{3}{2};
OCMATAE.dimplicitcontroldx=funch{4};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
OCMATAE.bctransversality=funch{5}{3};
OCMATAE.equilibrium=funch{5}{4};
if ~isempty(userbc) && userbc
    OCMATAE.userfuncbc=funch{5}{7};
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
if OCMATAE.objectivevaluecalc
    OCMATAE.objectivefunction=funch{8}{1};
    OCMATAE.objectivefunctionjacobian=funch{8}{2};
    OCMATAE.objectivefunctionparameterjacobian=funch{8}{3};
    if ~isautonomous(ocObj)
        OCMATAE.objectivefunctionderivativetime=funch{8}{4};
    end
end
if ~OCMATAE.autonomous
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

if OCMATAE.hitfunction
    OCMATAE.targetfunction=funch{10}{1};
end
% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

%
% End function definition
%

% equilibrium information

ocEP=limitset(ocgAsym);
depvarEP=dependentvar(ocEP);
if isempty(asymptoticmatrix)
    JacobianMatrix=jacobian(ocEP,1);
    [OCMATAE.asymptoticmatrix,OCMATAE.numstable,OCMATAE.numunstable,OCMATAE.numcenter,infoStruct]=asymptoticbc(JacobianMatrix{1},OCMATAE.pathtype,'c',OCMATCONT.ZeroDeviationTolerance,OCMATCONT.AsymptoticBCMethod);
end
OCMATAE.EP.saddlepoint=depvarEP{1};
OCMATAE.EP.arcarg=arcargument(ocEP);
OCMATAE.EP.linearization=linearization(ocEP);
OCMATAE.EP.implicitcontrolnum=length(implicitcontrolcoordinate(ocObj,OCMATAE.EP.arcarg));
% generate and adapt solution structure sol

sol=odestruct(ocgAsym);

arcn=arcnum(ocgAsym);
numberofodes=odenumber(ocgAsym);

solveInfo=solverinfo(ocgAsym);
if OCMATAE.objectivevaluecalc && ~isfield(solveInfo,'objectivevaluecoordinate')
    arcpos=arcposition(ocgAsym);
    o=objectivefunction(ocObj,sol,1);
    o=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,sol,1)))];
    sol.y(end+1,:)=0;
    for ii=1:arcn
        sol.y(numberofodes(ii)+1,arcpos(1,ii):arcpos(2,ii))=o(arcpos(1,ii):arcpos(2,ii));
    end
    numberofodes=numberofodes+1;
    OCMATAE.objectivevaluecoordinate=numberofodes;
elseif isfield(solveInfo,'objectivevaluecoordinate') && OCMATAE.objectivevaluecalc 
    OCMATAE.objectivevaluecoordinate=solveInfo.objectivevaluecoordinate;
    
else
    OCMATAE.objectivevaluecoordinate=[];
end
if OCMATAE.exogenousfunction
    OCMATAE.exogenousdynamicscoordinate=zeros(OCMATAE.exogenousnumberofstates,arcn);
    for ii=1:arcn
        sol.y(numberofodes(ii)+(1:OCMATAE.exogenousnumberofstates),arcpos(1,ii):arcpos(2,ii))=zeros(OCMATAE.exogenousnumberofstates,arcpos(2,ii)-arcpos(1,ii)+1);
        OCMATAE.exogenousdynamicscoordinate(1:OCMATAE.exogenousnumberofstates,ii)=numberofodes+(1:numberofodes+OCMATAE.exogenousnumberofstates);
    end
else
    OCMATAE.exogenousdynamicscoordinate=[];
end

depvar=dependentvar(ocgAsym);
OCMATAE.initialpoint=depvar{1}(:,1);

OCMATAE.arcarg=arcargument(ocgAsym);
OCMATAE.edge=[OCMATAE.arcarg(1:end-1);OCMATAE.arcarg(2:end)];

OCMATAE.numarc=arcnum(ocgAsym);
OCMATAE.arccoordinate=1:OCMATAE.numarc;
OCMATAE.trajectoryindex=OCMATAE.arccoordinate;

if OCMATAE.fixdistance
    OCMATAE.distance=norm(OCMATAE.EP.saddlepoint(OCMATAE.fixdistancecoordinate)-depvar{end}(OCMATAE.fixdistancecoordinate,end));
end

OCMATAE.continuationvector=[];
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

        arcintv=arcinterval(ocgAsym);
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

        OCMATAE.depvarcoordinate=continuationcoordinate;
        OCMATAE.startvalue=depvar{1}(continuationcoordinate,1);
        if length(OCMATAE.targetvalue)~=length(OCMATAE.startvalue)
            ocmatmsg('Length of target value and continuation coordinate are different.\n')
            sol=[];
            return
        end
        OCMATAE.continuationvector=OCMATAE.targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target point


        % free

        continuationparameter=0;

    case 'parameter'
        % target valur for parameters is a row vector
        OCMATAE.targetvalue=OCMATAE.targetvalue(:).';
        if isempty(continuationparameter)
            ocmatmessage('No continuation parameter provided.\n')
            return
        end
        if isempty(fixinitstate)
            fixinitstate=statecoord(ocObj);
        end

        arcintv=arcinterval(ocgAsym);
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

        if ~OCMATAE.simple
            switch OCMATAE.pathtype
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
        OCMATAE.modelparameter=modelparameter(ocgAsym);
        if ~isempty(OCMATAE.targetvalue)
            OCMATAE.startvalue=OCMATAE.modelparameter(OCMATAE.continuationparameterindex);
            OCMATAE.continuationvector=OCMATAE.targetvalue-OCMATAE.startvalue; % continue solution along the line from starting to target value
            continuationparameter=0;
        else
            continuationparameter=OCMATAE.modelparameter(OCMATAE.continuationparameterindex);
        end
        if isempty(fixinitstate)
            OCMATAE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATAE.initialstatecoordinate=fixinitstate;
        end
    case 'time'
        if isempty(fixinitstate)
            OCMATAE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATAE.initialstatecoordinate=fixinitstate;
        end
        if isempty(targettype)
            targettype='T';
        end
        OCMATAE.targettype=targettype;

        arcintv=arcinterval(ocgAsym);
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
        % if OCMATAE.freeendtime
        %     sol.parameters=[sol.parameters OCMATAE.endtime];
        %     OCMATAE.endtimecoordinate=parametercoordinate+1;
        %     parametercoordinate=length(sol.parameters);
        % end

        % free
        if continuationtime==length(arcintv)-1
            continuationparameter=OCMATAE.endtime;
        else
            continuationparameter=OCMATAE.OCMATAE.switchtimes(continuationtime);
        end
        OCMATAE.continuationtimeindex=continuationtime;
        OCMATAE.continuationcoordinate=parametercoordinate+1;
end

OCMATAE.fixinitstate=fixinitstate;
OCMATAE.freeparameter=freeparameter;
if ~isempty(freeparameter)
    freeparidx=parameterindex(ocObj,freeparameter);
    OCMATAE.freeparameterindex=freeparidx;
    parametercoordinate=length(sol.parameters);
    sol.parameters=[sol.parameters parametervalue(ocObj,freeparidx)];
    OCMATAE.freeparametercoordinate=parametercoordinate+(1:length(OCMATAE.freeparameterindex));
    parametercoordinate=length(sol.parameters);

    sol.parameters=[sol.parameters OCMATAE.EP.saddlepoint.'];
    OCMATAE.EP.equilibriumcoordinate=parametercoordinate+(1:length(OCMATAE.EP.saddlepoint));
    parametercoordinate=length(sol.parameters);

    switch OCMATAE.pathtype
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

end
sol.parameters=[sol.parameters continuationparameter];
OCMATAE.continuationcoordinate=parametercoordinate+1;
if strcmp(continuationtype,'time') && OCMATAE.freeendtime
    OCMATAE.endtimecoordinate=parametercoordinate+1;
end


pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);


OCMATCONT.codimension=1;
OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);