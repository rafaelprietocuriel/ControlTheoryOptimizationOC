function sol=initocmat_FTE4implicit(ocObj,ocgTrj,continuationtype,varargin)
%
% initocmat_FTE4implicit initialization of a finite time path calculation
% with a model exhibiting implicit controls.
%
% SOL=INITOCMAT_FTE4IMPLICIT(OCOBJ,OCGTRJ,CONTINUATIONTYPE,...)
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
initialstate=[];
initialcostate=[];
initialendtime=[];
initialarcargument=[];
findoptimalhorizon=[];
hitstatevalue=[];
hitstatecoordinate=[];
fixendstate=[];
fixendcostate=[];
fixinitstate=[];
fixinitcostate=[];
objectivevaluecalc=[];
exogenousfunction=[];
jumpargument=[];
entryindex=[];
jumpid=[];
opt=[];
freeparameter=[];
fixtimeindex=[];
optimalhorizon=false;
userbc=[];
hitvalue=[];
timeindex=[];
targetvalue=[];
initialcoordinate=[];
continuationcoordinate=[];
continuationparameter=[];
continuationtime=[];
initialcontrol=[];
method='';
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end

for ii=1:2:nargin-3
    if ischar(varargin{ii})
        if strcmp(varargin{ii},'option')
            eval('opt=varargin{ii+1};')
        else
            eval([varargin{ii} '=varargin{ii+1};'])
        end
    end
end
if isempty(findoptimalhorizon)
    findoptimalhorizon=false;
end
if isempty(hitstatevalue)
    hitstatevalue=[];
    hitstatecoordinate=[];
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(method)
    method='zero';
end
OCMATFTE.exogenousfunction=exogenousfunction;
OCMATFTE.objectivevaluecalc=objectivevaluecalc;

statecostatecoordinate=1:2*statenum(ocObj);
OCMATFTE.statecostatecoordinate=statecostatecoordinate;

%%% generate initial octrajectory if not provided
if isempty(ocgTrj)
    if isempty(initialstate)
        ocmatmsg('Initial state value is missing.\n')
        return
    end
    initialstate=initialstate(:);
    if isempty(initialendtime)
        ocmatmsg('No initial end time provided. Assuming T=0 at the initial continuation step.\n')
        initialendtime=0;
    end
    if isempty(initialarcargument)
        ocmatmsg('No initial arcid provided. Assuming arcid=0.\n')
        initialarcargument=0;
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');
    n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocgTrj.x=linspace(0,1,n);
    initocPt.y=initialstate;
    initocPt.modelparameter=parametervalue(ocObj);
    initocPt.modelname=modelname(ocObj);
    initocPt.x=0;
    initocPt.arcarg=initialarcargument;
    initocPt.arcinterval=[0 initialendtime];
    initocPt.arcposition=[1;1];
    if isempty(initialcostate)
        initialcostate=transversalitycondition(ocgtrajectory(initocPt));
    end
    ocgTrj.y=[initialstate(:);initialcostate];
    coord=implicitcontrolcoordinate(ocObj,initialarcargument);
    if ~isempty(coord)
        if isempty(initialcontrol)
            ocmatmsg('Initial control value is missing.\n')
            return
        end
        initialcontrol=initialcontrol(:);
        initocPt=calcic(ocObj,[ocgTrj.y;initialcontrol],initialarcargument,opt,'method',method);
        infoS=solverinfo(initocPt{1});

        if MessageDisplay
            ocmatmsg('\nInfo returned by the equation solver.\n\n')
            disp(infoS)
        end
        if infoS.exitflag==1
            if MessageDisplay
                ocmatmsg('\nTest admissibility.\n\n')
            end
            [b,infoS]=testadmissibility(initocPt{1},ocObj,opt);
            if b
                for ii=1:length(infoS)
                    if MessageDisplay
                        ocmatmsg('\nReturned violation values.\n\n%s\n\n',num2str(infoS(ii).violationmat.'))
                        ocmatmsg('\nReturned constraint values.\n\n%s\n\n',num2str(infoS(ii).constraintvalue.',2))
                    end
                end
                ocmatmsg('Detected solution is not admissible.\n')
                return
            end
            ctrl=control(initocPt{1});
            ocgTrj.y=[ocgTrj.y;ctrl{1}(coord)];
        else
            ocmatmsg('Cannot find initial control value.\n')
            return
        end
    end
    ocgTrj.y=ocgTrj.y(:,ones(1,n));
    ocgTrj.arcposition=[1;n];
    ocgTrj.arcinterval=[0 initialendtime];
    ocgTrj.arcarg=initialarcargument;
    ocgTrj.x0=0;
    ocgTrj.timehorizon=initialendtime;
    ocgTrj.modelparameter=parametervalue(ocObj);
    ocgTrj.modelname=modelname(ocObj);
    ocgTrj=ocgtrajectory(ocgTrj);
end

if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousinitialstates=exogenousinitialstates;
    OCMATFTE.exogenousnumberofstates=length(exogenousinitialstates);
end


% test if pure state constraints are defined
OCMATFTE.stateconstraint=stateconstraint(ocObj)&&~isempty(jumpargument);
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4FiniteHorizonPathContinuation');

% initialize global variable (OCMATFTE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATFTE.canonicalsystem=funch{1};
OCMATFTE.canonicalsystemjacobian=funch{2}{1};
OCMATFTE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATFTE.canonicalsystemhessian=funch{3}{1};
OCMATFTE.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATFTE.bcinitial=funch{5}{1};
OCMATFTE.bctransversality=funch{5}{2};
OCMATFTE.bcoptimalhorizon=funch{5}{3};
if OCMATFTE.stateconstraint
    OCMATFTE.bcstateconstraint=funch{5}{7};
    OCMATFTE.bctransversalitysc=funch{5}{8};
end


% function for Jacobian
OCMATFTE.bcjacobianinitial=funch{6}{1};
OCMATFTE.bcjacobiantransversality=funch{6}{2};

% function describing the hybrid structure of the problem
OCMATFTE.hybridinfo=funch{7}{1};
OCMATFTE.domain=funch{7}{2};
OCMATFTE.guard=funch{7}{3};
OCMATFTE.reset=funch{7}{4};
OCMATFTE.switchtime=funch{7}{5};
OCMATFTE.jacobianguard=funch{7}{7};
OCMATFTE.jacobianreset=funch{7}{8};
OCMATFTE.domaindiscretization=funch{7}{9};

OCMATFTE.findoptimalhorizon=findoptimalhorizon;

if OCMATFTE.objectivevaluecalc
    OCMATFTE.objectivefunction=funch{8}{1};
    OCMATFTE.objectivefunctionjacobian=funch{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATFTE.objectivefunctionderivativetime=funch{8}{4};
    OCMATFTE.salvagevalue=funch{5}{6};
end
if OCMATFTE.exogenousfunction
    OCMATFTE.exogenousdynamics=funch{4}{1};
    OCMATFTE.exogenousjacobian=funch{4}{2};
    OCMATFTE.exogenousparameterjacobian=funch{4}{3};
end

% general function
OCMATFTE.plotcontinuation=funch{11};
OCMATFTE.testadmissibility=funch{12};
OCMATFTE.datapath=funch{20};
OCMATFTE.saveintermediatefiles=funch{21};

OCMATFTE.parametervalue=parametervalue(ocObj);
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATFTE.statecostatecoord=[scoord(:).' cscoord(:).'];
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
OCMATFTE.freeparameter=freeparameter;
if ~isempty(freeparameter)
    OCMATFTE.freeparameter=parameterindex(ocObj,freeparameter);
end
% mode and path specific variables
depvar=dependentvar(ocgTrj);
OCMATFTE.parametervalue=parametervalue(ocObj);
OCMATFTE.autonomous=isautonomous(ocObj);


OCMATFTE.initialcoordinate=initialcoordinate;
OCMATFTE.varyparameterindex=timeindex;
OCMATFTE.targetvalue=targetvalue;
OCMATFTE.optimalhorizon=optimalhorizon;

sol=odestruct(ocgTrj);
OCMATFTE.initialtime=sol.x0;

arcn=arcnum(ocgTrj);
numberofodes=zeros(1,length(depvar));
for ii=1:arcn
    numberofodes(ii)=size(depvar{ii},1);
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

        arcintv=arcinterval(ocgTrj);
        OCMATFTE.switchtimes=arcintv(2:end-1);
        OCMATFTE.endtime=arcintv(end);

        OCMATFTE.fixtimeindex=fixtimeindex;
        OCMATFTE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATFTE.switchtimes(OCMATFTE.freetimeindex);
        OCMATFTE.freetimecoordinate=parametercoordinate+(1:length(OCMATFTE.freetimeindex));
        OCMATFTE.switchtimecoordinate(OCMATFTE.freetimeindex)=OCMATFTE.freetimecoordinate;
        parametercoordinate=length(sol.parameters);

        OCMATFTE.depvarcoordinate=continuationcoordinate;
        OCMATFTE.startvalue=depvar{1}(continuationcoordinate,1);
        if length(OCMATFTE.targetvalue)~=length(OCMATFTE.startvalue)
            ocmatmsg('Length of target value and continuation coordinate are different.\n')
            sol=[];
            return
        end
        OCMATFTE.continuationvector=OCMATFTE.targetvalue-OCMATFTE.startvalue; % continue solution along the line from starting to target point


        % free

        continuationparameter=0;
    case 'parameter'
        OCMATFTE.targetvalue=OCMATFTE.targetvalue(:).';
        if isempty(continuationparameter)
            ocmatmessage('No continuation parameter provided.\n')
            return
        end
        if isempty(fixinitstate)
            fixinitstate=statecoord(ocObj);
        end

        arcintv=arcinterval(ocgTrj);
        OCMATFTE.switchtimes=arcintv(2:end-1);
        OCMATFTE.endtime=arcintv(end);

        OCMATFTE.fixtimeindex=fixtimeindex;
        OCMATFTE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATFTE.switchtimes(OCMATFTE.freetimeindex);
        OCMATFTE.freetimecoordinate=parametercoordinate+(1:length(OCMATFTE.freetimeindex));
        OCMATFTE.switchtimecoordinate(OCMATFTE.freetimeindex)=OCMATFTE.freetimecoordinate;
        parametercoordinate=length(sol.parameters);


        OCMATFTE.continuationparameterindex=parameterindex(ocObj,continuationparameter);
        OCMATFTE.freeparameterindex=parameterindex(ocObj,continuationparameter);
        OCMATFTE.modelparameter=modelparameter(ocgTrj);
        if ~isempty(OCMATFTE.targetvalue)
            OCMATFTE.startvalue=OCMATFTE.modelparameter(OCMATFTE.continuationparameterindex);
            OCMATFTE.continuationvector=OCMATFTE.targetvalue-OCMATFTE.startvalue; % continue solution along the line from starting to target value
            continuationparameter=0;
        else
            continuationparameter=OCMATFTE.modelparameter(OCMATFTE.continuationparameterindex);
        end
        if isempty(fixinitstate)
            OCMATFTE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATFTE.initialstatecoordinate=fixinitstate;
        end
    case 'time'
        if isempty(fixinitstate)
            OCMATFTE.initialstatecoordinate=statecoord(ocObj);
        else
            OCMATFTE.initialstatecoordinate=fixinitstate;
        end

        arcintv=arcinterval(ocgTrj);
        OCMATFTE.switchtimes=arcintv(2:end-1);
        OCMATFTE.endtime=arcintv(end);

        OCMATFTE.fixtimeindex=fixtimeindex;

        if isempty(continuationtime)
            continuationtime=length(arcintv)-1;
        end
        if continuationtime==length(arcintv)-1
            OCMATFTE.freetimeindex=setdiff(1:length(arcintv)-2,fixtimeindex);

        else
            OCMATFTE.freetimeindex=setdiff(1:length(arcintv)-2,[fixtimeindex continuationtime]);
        end


        parametercoordinate=length(sol.parameters);
        sol.parameters=OCMATFTE.switchtimes(OCMATFTE.freetimeindex);
        OCMATFTE.freetimecoordinate=parametercoordinate+(1:length(OCMATFTE.freetimeindex));
        OCMATFTE.switchtimecoordinate(OCMATFTE.freetimeindex)=OCMATFTE.freetimecoordinate; % coordinate of the switching times within the free parameters for each arc, 0 ... switching time is fixed
        parametercoordinate=length(sol.parameters);
        % if OCMATFTE.freeendtime
        %     sol.parameters=[sol.parameters OCMATFTE.endtime];
        %     OCMATFTE.endtimecoordinate=parametercoordinate+1;
        %     parametercoordinate=length(sol.parameters);
        % end

        % free
        if continuationtime==length(arcintv)-1
            continuationparameter=OCMATFTE.endtime;
        else
            continuationparameter=OCMATFTE.OCMATFTE.switchtimes(continuationtime);
        end
        OCMATFTE.continuationtimeindex=continuationtime;
        OCMATFTE.continuationcoordinate=parametercoordinate+1;
    otherwise
end
OCMATFTE.initialstate=depvar{1}(OCMATFTE.initialstatecoordinate,1);


OCMATFTE.transversalityconditioncs=0;
OCMATFTE.fixinitstate=fixinitstate;
OCMATFTE.freeparameter=freeparameter;
if ~isempty(freeparameter)
    freeparidx=parameterindex(ocObj,freeparameter);
    OCMATFTE.freeparameterindex=freeparidx;
    parametercoordinate=length(sol.parameters);
    sol.parameters=[sol.parameters parametervalue(ocObj,freeparidx)];
    OCMATFTE.freeparametercoordinate=parametercoordinate+(1:length(OCMATFTE.freeparameterindex));
    parametercoordinate=length(sol.parameters);
end


sol.parameters=[sol.parameters continuationparameter];
OCMATFTE.continuationcoordinate=parametercoordinate+1;

 

OCMATFTE.statenum=statenum(ocObj);
OCMATFTE.fixinitstatecoord=fixinitstate;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitcostatecoord=fixinitcostate;
OCMATFTE.fixendcostatecoord=fixendcostate;
OCMATFTE.initstate=depvar{1}(fixinitstate,1);
OCMATFTE.endstate=depvar{1}(fixendstate,end);
OCMATFTE.endcostate=depvar{1}(OCMATFTE.statenum+fixendstate,end);
OCMATFTE.hitstatevalue=hitstatevalue;
OCMATFTE.hitstatecoordinate=hitstatecoordinate;

OCMATFTE.variationalcalculation=[];
%OCMATFTE.continitstate=false;
% if isempty(fixinitstate)
%     OCMATFTE.startvalue=depvar(fixinitstate,1);
%     %OCMATFTE.continitstate=true;
% elseif isempty(fixendstateidx)
%     OCMATFTE.startvalue=depvar(fixendstate,end);
% end
% if ~isempty(fixendstate)
%     OCMATFTE.endstate=depvar(fixendstate,end);
%     OCMATFTE.initstate=[];
% elseif ~isempty(fixinitstate)
%     OCMATFTE.initstate=depvar(fixinitstate,1);
%     OCMATFTE.endstate=[];
% end
pathname=OCMATFTE.datapath();
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(initialcoordinate);

OCMATCONT.codimension=1;
OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);
