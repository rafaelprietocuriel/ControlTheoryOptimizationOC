function sol=initocdae_FTE(ocObj,ocTrj,conttype,continuationindex,varargin)
clear global OCMATFTE OCMATCONT 

global OCMATFTE OCMATCONT
sol=[];
initialstate=[];
initialcostate=[];
fixinitialstate=[];
fixendstate=[];
initialarcargument=[];
initialtime=[];
initialendtime=[];
targetvalue=[];
opt=[];
initcontpar=[];
collocationnum=[];
collocationmeth=[];
interpolationmeth=[];
objectivevaluecalc=[];
freeparameter=[];
timegrid=[];
recalculatecoeff=0;

initialstateidx=find(strcmpi(varargin,'initialstate'));
initialcostateidx=find(strcmpi(varargin,'initialcostate'));
initialendtimeidx=find(strcmpi(varargin,'initialendtime'));
initialtimeidx=find(strcmpi(varargin,'initialtime'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
optionidx=find(strcmpi(varargin,'option'));
initcontparidx=find(strcmpi(varargin,'initcontpar'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
collocationnumidx=find(strcmpi(varargin,'collocationnum'));
collocationmethidx=find(strcmpi(varargin,'collocationmethod'));
interpolationmethidx=find(strcmpi(varargin,'interpolationmethod'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
timegrididx=find(strcmpi(varargin,'timegrid'));
recalculatecoeffidx=find(strcmpi(varargin,'recalculatecoeff'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(initcontparidx)
    initcontpar=varargin{initcontparidx+1}; % only used if ocEP is empty
end
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if ~isempty(collocationnumidx)
    collocationnum=varargin{collocationnumidx+1};
end
if ~isempty(collocationmethidx)
    collocationmeth=varargin{collocationmethidx+1};
end
if ~isempty(interpolationmethidx)
    interpolationmeth=varargin{interpolationmethidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(timegrididx)
    timegrid=varargin{timegrididx+1};
end
if ~isempty(recalculatecoeffidx)
    recalculatecoeff=varargin{recalculatecoeffidx+1};
end
if isempty(initcontpar)
    initcontpar=0;
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(collocationmeth) && isempty(ocTrj)
    collocationmeth=getocoptions(opt,'SBVPOC','CollocationMethod');
elseif isempty(collocationmeth)
    collocationmeth=collocationmethod(ocTrj);
end
if isempty(collocationnum) && isempty(ocTrj)
    collocationnum=getocoptions(opt,'SBVPOC','CollocationNumber');
elseif isempty(collocationmeth)
    collocationnum=collocationnumber(ocTrj);
end
if isempty(interpolationmeth) && isempty(ocTrj)
    interpolationmeth=getocoptions(opt,'SBVPOC','InterpolationMethod');
elseif isempty(collocationmeth)
    interpolationmeth=interpolationmethod(ocTrj);
end
if ischar(freeparameter)
    freeparameter=parameterindex(ocObj,freeparameter);
end
% get modelspecific functions
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4DAEContinuation');
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATFTE.canonicalsystem=funch{1}{1};
OCMATFTE.dlagrangefunctiondu=funch{1}{2};
OCMATFTE.complementaryslacknesscondition=funch{1}{3};

OCMATFTE.canonicalsystemjacobian=funch{2}{1};
OCMATFTE.canonicalsystemparameterjacobian=funch{2}{2};
OCMATFTE.canonicalsystemhessian=funch{2}{3};
OCMATFTE.canonicalsystemparameterhessian=funch{2}{4};
OCMATFTE.dlagrangefunctiondujacobian=funch{2}{5};
OCMATFTE.dlagrangefunctionduparameterjacobian=funch{2}{6};
OCMATFTE.complementaryslacknessconditionjacobian=funch{2}{7};
OCMATFTE.complementaryslacknessconditionparameterjacobian=funch{2}{8};

OCMATFTE.constraint=funch{3}{1};

OCMATFTE.bcinitial=funch{5}{1};
OCMATFTE.bctransversality=funch{5}{2};
OCMATFTE.salvagevalue=funch{5}{3};

OCMATFTE.objectivefunction=funch{8}{1};
OCMATFTE.objectivefunctionjacobian=funch{8}{2};
OCMATFTE.objectivefunctionparameterjacobian=funch{8}{3};
%OCMATFTE.objectivefunctionderivativetime=funch{8}{4};

OCMATFTE.plotcontinuation=funch{11};
OCMATFTE.datapath=funch{20};
OCMATFTE.saveintermediatefiles=funch{21};

%%% generate initial octrajectory if not provided
if isempty(ocTrj)
    if ~isempty(initialstateidx)
        initialstate=varargin{initialstateidx+1};
        initialstate=initialstate(:);
    end
    if ~isempty(initialcostateidx)
        initialcostate=varargin{initialcostateidx+1};
        initialcostate=initialcostate(:);
    end
    if ~isempty(initialendtimeidx)
        initialendtime=varargin{initialendtimeidx+1};
    end
    if ~isempty(initialtimeidx)
        initialtime=varargin{initialtimeidx+1};
    end
    if ~isempty(initialarcargumentidx)
        initialarcargument=varargin{initialarcargumentidx+1};
    end
    if ~isempty(optionidx)
        opt=varargin{optionidx+1};
    end
    if isempty(opt)
        opt=defaultocoptions;
    end
    if isempty(initialtime)
        initialtime=0;
    end
    if isempty(initialendtime)
        initialendtime=0;
    end
    if isempty(initialarcargument)
        initialarcargument=0;
    end
    initocPt.y=initialstate;
    initocPt.x=0;
    initocPt.arcarg=initialarcargument;
    initocPt.arcinterval=[0 initialendtime];
    initocPt.arcposition=[1;1];
    if isempty(initialcostate)
        initialcostate=transversalitycondition(ocObj,octrajectory(initocPt));
    end
    initocPt.y=[initialstate;initialcostate];
    initialcontrol=control(ocObj,octrajectory(initocPt));
    initiallm=lagrangemultiplier(ocObj,octrajectory(initocPt));

    % create initial octrajectory
    n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocTrj.x=linspace(0,1,n);
    ocTrj.y=repmat([initialstate;initialcostate;initialcontrol;initiallm],1,n);
    ocTrj.t0=initialtime;
    ocTrj.arcarg=initialarcargument;
    ocTrj.arcinterval=[0 1];
    ocTrj.timehorizon=initialendtime;
    ocTrj.modelparameter=parametervalue(ocObj);
    ocTrj.modelname=modelname(ocObj);
    ocTrj.solverinfo.type='dae';
    ocTrj.solverinfo.daeorder=[ones(2*statenum(ocObj),1);zeros(controlnum(ocObj),1);zeros(controlconstraintnum(ocObj),1)];
    ocTrj.solverinfo.statecoordinate=1:statenum(ocObj);
    ocTrj.solverinfo.costatecoordinate=statenum(ocObj)+(1:statenum(ocObj));
    ocTrj.solverinfo.controlcoordinate=ocTrj.solverinfo.costatecoordinate(end)+(1:controlnum(ocObj));
    ocTrj.solverinfo.lagrangemultipliercoordinate=ocTrj.solverinfo.controlcoordinate(end)+(1:controlconstraintnum(ocObj));
    ocTrj=octrajectory(ocTrj);
    ocTrj=octrajectory2ocdae(ocTrj);
end

timemesh=independentvar(ocTrj);
ordr=daeorder(ocTrj);
objectivevaluecoord=objectivevaluecoordinate(ocTrj);
if objectivevaluecalc
    if isempty(objectivevaluecoord)
        ordr=[ordr;1];
        ocTrj.solverinfo.daeorder=ordr;
        ocTrj.solverinfo.objectivevaluecoordinate=length(ordr);
        OT=discountedsalvagevalue(ocObj,ocTrj);
        o=objectivefunction(ocObj,ocTrj,1);
        ocTrj.y(length(ordr),:)=OT+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(timemesh))];
        objectivevaluecoord=length(ordr);
    end
else
    if ~isempty(objectivevaluecoord)
        ordr(objectivevaluecoord)=[];
        ocTrj.y(objectivevaluecoord,:)=[];
    end
end
OCMATFTE.objectivevaluecalc=objectivevaluecalc;

sol.x=timemesh;
sol.y=dependentvar(ocTrj);
if ~recalculatecoeff
    sol.data.coeff=rungekuttacoefficient(ocTrj);
end
sol.data.xcol=collocationmesh(ocTrj);

OCMATFTE.statenum=statenum(ocObj);

OCMATFTE.statecoordinate=statecoordinate(ocTrj);
OCMATFTE.costatecoordinate=costatecoordinate(ocTrj);
OCMATFTE.controlcoordinate=controlcoordinate(ocTrj);
OCMATFTE.lagrangemultipliercoordinate=lagrangemultipliercoordinate(ocTrj);
OCMATFTE.totalcoordinate=1:(2*statenum(ocObj)+controlnum(ocObj)+inequalitycontrolconstraintnum(ocObj));
OCMATFTE.objectivevaluecoordinate=objectivevaluecoord;

if isempty(fixinitialstate) % default: initial states are fixed and end states are free
    fixinitialstate=1:OCMATFTE.statenum;
end
freeinitialstate=setdiff(1:OCMATFTE.statenum,fixinitialstate);
freeendstate=setdiff(1:OCMATFTE.statenum,fixendstate);

% filenames for storing intermediate continuation results and possible
% animation file
pathname=OCMATFTE.datapath();
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

%
OCMATFTE.basicmoviefilename=fullocmatfile(pathname,'AnimCont');

OCMATFTE.modelname=modelname(ocObj);
OCMATFTE.modelparameter=parametervalue(ocObj);
OCMATFTE.freeparametercoordinate=[];
OCMATFTE.freeparameterindex=freeparameter;
OCMATFTE.initialstate=state(ocObj,ocTrj);
OCMATFTE.initialstate=OCMATFTE.initialstate(:,1);
if ~isempty(modelname(ocTrj)) && ~strcmp(modelname(ocObj),modelname(ocTrj))
    ocmatmsg('Warning, model name of extremal and model are different')
end

OCMATFTE.initialarctimes=arcinterval(ocTrj);
% determine the free parameters
sol.parameters=[];
if ~isempty(freeparameter)
    sol.parameters=[sol.parameters parametervalue(ocObj,freeparameter)];
    OCMATFTE.freeparametercoordinate=1:length(freeparameter);
end
freeparameterctr=length(sol.parameters);
% continuation type:  time,
OCMATFTE.conttype=conttype;
switch conttype
    case 'time'
        if isempty(continuationindex)
            OCMATFTE.continuationindex=2; % continue endtime
        else
            OCMATFTE.continuationindex=continuationindex; % 1 initial time, 2 endtime
        end
    case 'initialstate'
        if isempty(targetvalue)
            return
        end
        if length(targetvalue)~=length(continuationindex)
            return
        end
    case 'endstate'
        if isempty(targetvalue)
            return
        end
        if length(targetvalue)~=length(continuationindex)
            return
        end
    case 'parameter'
        if ischar(continuationindex)
            continuationindex=parameterindex(ocObj,continuationindex);
        end
        OCMATFTE.continuationindex=continuationindex;
    otherwise
        return
end
OCMATFTE.targetvalue=targetvalue(:);
OCMATCONT.TargetValueNum=length(OCMATFTE.targetvalue);

switch conttype
    case 'time'
        if OCMATFTE.continuationindex==1
            %sol.parameters=[sol.parameters initialtime(ocTrj)];
            continuationparametervalue=initialtime(ocTrj)+initcontpar;
        elseif OCMATFTE.continuationindex==2
            %sol.parameters=[sol.parameters timehorizon(ocTrj)+initcontpar];
            sol.x=(sol.x-sol.x(1))/sol.x(end);
            continuationparametervalue=timehorizon(ocTrj)+initcontpar;
        else
            return
        end
    case 'initialstate'
        %sol.parameters=0;
        continuationparametervalue=initcontpar;
        fixinitialstate=continuationindex;
        OCMATFTE.initialstatecoordinate=continuationindex;
        OCMATFTE.targetvector=OCMATFTE.targetvalue-OCMATFTE.initialstate(continuationindex);
    case 'endstate'
        %sol.parameters=0;
        continuationparametervalue=initcontpar;
        fixendstate=continuationindex;
        OCMATFTE.endstate=sol.y(OCMATFTE.statecoordinate,end);

        OCMATFTE.endstatecoordinate=continuationindex;
        OCMATFTE.targetvector=OCMATFTE.targetvalue-OCMATFTE.endstate(continuationindex);
    case 'parameter'
        continuationparametervalue=parametervalue(ocObj,continuationindex); 
        if ~isempty(OCMATFTE.targetvalue)
            %    OCMATFTE.targetvector=OCMATFTE.targetvalue-OCMATFTE.startingparametervalue;
        else
            %    OCMATFTE.targetvector=1;
        end
end

% coordinate is used as a signifier for the position (row) in the solution
% path
% index is uses as a signifier for the position in the free parameters
% vector
OCMATFTE.fixinitialstatecoordinate=fixinitialstate;
OCMATFTE.freeinitialstatecoordinate=freeinitialstate;
OCMATFTE.fixendstatecoordinate=fixendstate;
OCMATFTE.freeendstatecoordinate=freeendstate;

OCMATFTE.fixinitialstate=sol.y(OCMATFTE.fixinitialstatecoordinate,1);
OCMATFTE.fixendstate=sol.y(OCMATFTE.fixendstatecoordinate,end);

if ~isempty(fixendstate)
    OCMATFTE.freeendcostatecoordinate=fixendstate;
    freeparameterctr=freeparameterctr+(1:length(fixendstate));
    sol.parameters=[sol.parameters sol.y(OCMATFTE.costatecoordinate(fixendstate),end).'];
    OCMATFTE.freeendcostateindex=fixendstate;
end

if ~isempty(freeinitialstate)
    OCMATFTE.fixinitialcostatecoordinate=freeinitialstate;
    freeparameterctr=freeparameterctr+(1:length(fixendstate));
    sol.parameters=[sol.parameters sol.y(OCMATFTE.costatecoordinate(freeinitialstate),1).'];
    OCMATFTE.fixinitialcostateindex=freeinitialstate;
end
sol.parameters=[sol.parameters continuationparametervalue];
if ~isempty(timegrid)
    sol.xnew=timegrid;
end
%tic
sol=adapt2mdaesolver(sol,ordr,collocationmeth,collocationnum,interpolationmeth);
%t=toc;fprintf('elapsed time: %f\n',t)

OCMATFTE.continuationcoordinate=freeparameterctr(end)+1;
OCMATFTE.continuationtype=conttype;

% Structure Jacobian
OCMATFTE.dFDPAR=zeros(OCMATCONT.componentnumber,length(sol.parameters));
OCMATFTE.dFDX=zeros(OCMATCONT.componentnumber);

