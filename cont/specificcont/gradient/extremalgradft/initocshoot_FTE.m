function sol=initocshoot_FTE(ocObj,ocgTrj,conttype,contindex,varargin)
clear global OCGRADSOL OCSHOOTCONT
global OCGRADSOL OCSHOOTCONT
sol=[];
initstate=[];
fixinitialstate=[];
fixendstate=[];
initialcostate=[];
initialarcargument=[];
initialtime=[];
initialendtime=[];
targetvalue=[];
objoctype='';
opt=[];
initcontpar=[];

initstateidx=find(strcmpi(varargin,'initialstate'));
initialcostateidx=find(strcmpi(varargin,'initialcostate'));
initialendtimeidx=find(strcmpi(varargin,'initialendtime'));
initialtimeidx=find(strcmpi(varargin,'initialtime'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
octypeidx=find(strcmpi(varargin,'octype'));
optionidx=find(strcmpi(varargin,'option'));
initcontparidx=find(strcmpi(varargin,'initcontpar'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(initcontparidx)
    initcontpar=varargin{initcontparidx+1}; % only used if ocEP is empty
end
if ~isempty(octypeidx)
    objoctype=varargin{octypeidx+1}; % only used if ocgTrj is empty
end
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if isempty(initcontpar)
    initcontpar=0;
end
% get modelspecific functions
OCSHOOTCONT.modelfunc=modelspecificfunc(ocObj,'4ShootFiniteHorizonPathContinuation');
funch=OCSHOOTCONT.modelfunc(); % model specific function handles for saddle path continuation

OCGRADSOL.canonicalsystemdynamics=funch{1};
OCGRADSOL.objectivefunction=funch{2};
OCGRADSOL.disountedsalvagevalue=funch{3};
OCGRADSOL.transversaltycondition=funch{4};
OCGRADSOL.extendedhamiltonian=funch{5};
OCGRADSOL.lagrangemultiplier=funch{6};
OCGRADSOL.admissible=funch{7};
OCGRADSOL.optimalcontrol=funch{8};
OCGRADSOL.constraint=funch{9};

OCGRADSOL.plotcontinuation=funch{11};
OCGRADSOL.datapath=funch{20};
OCGRADSOL.saveintermediatefiles=funch{21};

%%% generate initial octrajectory if not provided
if isempty(ocgTrj)
    if ~isempty(initstateidx)
        initstate=varargin{initstateidx+1};
        initstate=initstate(:);
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
    if isempty(objoctype)
        objoctype='concentrated'; % default is local (not distributed)
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
    initocPt.y=initstate;
    initocPt.x=0;
    initocPt.arcarg=initialarcargument;
    initocPt.arcinterval=[0 initialendtime];
    initocPt.arcposition=[1;1];
    if isempty(initialcostate)
        initialcostate=transversalitycondition(ocObj,octrajectory(initocPt));
    end
    initocPt.y=[initstate;initialcostate];
    ctrl=control(ocObj,octrajectory(initocPt));

    % create initial ocgradtrajectory
    n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocgTrj.argument.t=linspace(0,1,n);
    ocgTrj.variable.y=initstate(:);
    ocgTrj.variable.y=ocgTrj.variable.y(:,ones(1,n));
    ocgTrj.variable.cst_y=initialcostate;
    ocgTrj.variable.cst_y=ocgTrj.variable.cst_y(:,ones(1,n));
    ocgTrj.variable.v=ctrl;
    ocgTrj.variable.v=ocgTrj.variable.v(:,ones(1,n));
    ocgTrj.variable.o=0;
    ocgTrj.variable.o=ocgTrj.variable.o(:,ones(1,n));
    ocgTrj.variable.arcid=repmat(initialarcargument,1,n);
    ocgTrj.initargument.t0=initialtime;
    ocgTrj.timehorizon=initialendtime;
    ocgTrj.modelparameter=parametervalue(ocObj);
    ocgTrj.modelname=modelname(ocObj);
    ocgTrj.type=objoctype;
    ocgTrj.objectivevalue=OCGRADSOL.disountedsalvagevalue(initialendtime,initstate,ocgTrj.modelparameter);
    ocgTrj=ocgradtrajectory(ocgTrj);

end
cst_y=costate(ocgTrj);
y=state(ocgTrj);
% number of states and controls
[control_num,state_num]=numberofvariables(ocObj);

OCSHOOTCONT.statecoordinate=1:state_num.concentrated;
OCSHOOTCONT.costatecoordinate=state_num.concentrated+(1:state_num.concentrated);
OCSHOOTCONT.objectivecoordinate=OCSHOOTCONT.costatecoordinate(end)+1;
if isempty(fixinitialstate) % default: initial states are fixed and end states are free
    fixinitialstate=1:state_num.concentrated;
end
freeinitialstate=setdiff(1:state_num.concentrated,fixinitialstate);
freeendstate=setdiff(1:state_num.concentrated,fixendstate);

% coordinate is used as a signifier for the position (row) in the solution
% path
% index is uses as a signifier for the position in the free parameters
% vector
OCSHOOTCONT.fixinitialstatecoordinate=fixinitialstate;
OCSHOOTCONT.freeinitialstatecoordinate=freeinitialstate;
OCSHOOTCONT.fixendstatecoordinate=fixendstate;
OCSHOOTCONT.freeendstatecoordinate=freeendstate;

OCSHOOTCONT.fixinitialstate=y(OCSHOOTCONT.fixinitialstatecoordinate,1);
OCSHOOTCONT.fixendstate=y(OCSHOOTCONT.fixendstatecoordinate,end);

% information on control, mixed-control constraint
cstrnum=controlconstraintnum(ocObj);
mixedcstr=mixedconstraint(ocObj);

OCSHOOTCONT.cconstraint=0;
OCSHOOTCONT.mcconstraint=0;
if cstrnum
    OCSHOOTCONT.cconstraint=1;
    if mixedcstr
        OCSHOOTCONT.mcconstraint=1;
    else
        OCSHOOTCONT.mcconstraint=0;
    end
end
% filenames for storing intermediate continuation results and possible
% animation file
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

% 
OCSHOOTCONT.basicmoviefilename=fullocmatfile(pathname,'AnimCont');

OCSHOOTCONT.modelname=modelname(ocObj);
OCSHOOTCONT.modelparameter=parametervalue(ocObj);

if ~isempty(modelname(ocgTrj)) && ~strcmp(modelname(ocObj),modelname(ocgTrj))
    ocmatmsg('Warning, model name of extremal and model are different')
end

% determine the free parameters
freeparameterctr=0;
sol.parameters=y(OCSHOOTCONT.freeinitialstatecoordinate,1).'; % free initial states have to be determined by some equations
OCSHOOTCONT.freeinitialstateindex=(freeparameterctr+1):length(sol.parameters);
freeparameterctr=length(sol.parameters);
sol.parameters=[sol.parameters cst_y(OCSHOOTCONT.fixinitialstatecoordinate,1).'];
OCSHOOTCONT.fixinitialstateindex=(freeparameterctr+1):length(sol.parameters);
freeparameterctr=length(sol.parameters);

% continuation type:  time,
OCSHOOTCONT.conttype=conttype;
switch conttype
    case 'time'
        if isempty(contindex)
            OCSHOOTCONT.contindex=2; % continue endtime
        else
            OCSHOOTCONT.contindex=contindex; % 1 initial time, 2 endtime
        end
    case 'initialstate'
        if isempty(targetvalue)
            return
        end
        if length(targetvalue)~=length(contindex)
            return
        end
    case 'parameter'
        if ischar(contindex)
            contindex=parameterindex(ocObj,contindex);
        end
        OCSHOOTCONT.contindex=contindex;
    otherwise
        return
end

% OCSHOOTCONT.continuationindex=1;
% fn=fieldnames(control_num);
% laststatecoord=0;
% for ii=1:length(fn)
%     if ~isempty(state_num.(fn{ii}))
%         OCSHOOTCONT.statecoordinate.(fn{ii})=1:state_num.(fn{ii});
%         laststatecoord=laststatecoord+state_num.(fn{ii});
%     end
% end
% OCSHOOTCONT.objectivecoordinate=laststatecoord+1;
OCSHOOTCONT.control_num=control_num;
OCSHOOTCONT.state_num=state_num;
OCSHOOTCONT.targetvalue=targetvalue(:);
OCSHOOTCONT.octype=octype(ocgTrj);

sol.extremal=generategradextremal(ocgTrj);
switch conttype
    case 'time'
        if OCSHOOTCONT.contindex==1
            sol.parameters=[sol.parameters initialtime(ocgTrj)];
        elseif OCSHOOTCONT.contindex==2
            sol.parameters=[sol.parameters timehorizon(ocgTrj)+initcontpar];
            sol.extremal.t=(sol.extremal.t-sol.extremal.t0)/sol.extremal.t(end);
        else
            return
        end
    case 'initialstate'
        sol.parameters=[sol.parameters 0];
        OCSHOOTCONT.startingstatevalue=initialstate(ocgTrj);
        OCSHOOTCONT.targetvector=OCSHOOTCONT.targetvalue-OCSHOOTCONT.startingstatevalue;
    case 'parameter'
        sol.parameters=[sol.parameters parametervalue(ocObj,contindex)];
        %OCSHOOTCONT.startingparametervalue=parametervalue(ocObj,contindex);
        if ~isempty(OCSHOOTCONT.targetvalue)
        %    OCSHOOTCONT.targetvector=OCSHOOTCONT.targetvalue-OCSHOOTCONT.startingparametervalue;
        else
        %    OCSHOOTCONT.targetvector=1;
        end
end
OCSHOOTCONT.continuationindex=freeparameterctr+1;
switch OCSHOOTCONT.octype
    case 'concentrated'
        n=length(time(ocgTrj));
        OCSHOOTCONT.loccontrolcoordinate=reshape(1:control_num.concentrated*n,[],n);
        OCSHOOTCONT.initiallocstatecoordinate=OCSHOOTCONT.loccontrolcoordinate(end)+(1:state_num.concentrated).';
        OCSHOOTCONT.endloccostatecoordinate=OCSHOOTCONT.initiallocstatecoordinate(end)+(1:state_num.concentrated).';
end
OCSHOOTCONT.TIMEMESH.num=length(sol.extremal.t);

function [control_num,state_num]=numberofvariables(ocObj)

control_num.concentrated=[];
control_num.dist=[];
control_num.agg=[];

state_num.concentrated=[];
state_num.dist=[];
state_num.agg=[];
if isocmodel(ocObj)
    control_num.concentrated=controlnum(ocObj);
    state_num.concentrated=statenum(ocObj);
else
end

