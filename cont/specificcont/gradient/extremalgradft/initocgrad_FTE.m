function sol=initocgrad_FTE(ocObj,ocgTrj,conttype,contindex,varargin)
clear global OCGRADSOL OCGRADCONT
global OCGRADSOL OCGRADCONT
sol=[];
initstate=[];
initialcostate=[];
fixinitialstate=[];
fixendstate=[];
initialarcargument=[];
initialtime=[];
initialendtime=[];
targetvalue=[];
objoctype='';
opt=[];
initcontpar=[];
method='';

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
methodidx=find(strcmpi(varargin,'method'));
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
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if ~isempty(methodidx)
    method=varargin{methodidx+1};
end

if isempty(initcontpar)
    initcontpar=0;
end
if isempty(method)
    method='grad';
end

if ~strcmp(method,'grad')
    ocmatmsg('Only ''grad'' method is implemented yet')
end
% get modelspecific functions
OCGRADCONT.modelfunc=modelspecificfunc(ocObj,'4GradContinuation');
funch=OCGRADCONT.modelfunc(); % model specific function handles for saddle path continuation

OCGRADSOL.statedynamics=funch{1}{1};
OCGRADSOL.costatedynamics=funch{2}{1};
OCGRADSOL.objectivefunction=funch{3}{1};
OCGRADSOL.gradhamiltonian=funch{4}{1};
OCGRADSOL.salvagevalue=funch{5}{1};
OCGRADSOL.transversaltycondition=funch{6}{1};
OCGRADSOL.statedynamicsext=funch{7}{1};
OCGRADSOL.explicitgradientcontrol=funch{8}{1};
OCGRADSOL.projectionlocal=funch{9}{1};
OCGRADSOL.admissible=funch{10};
OCGRADSOL.constraint=funch{12};
OCGRADSOL.controlbounds=funch{13};
OCGRADSOL.globalcorrector=funch{16};

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
    ocgTrj.initargument.t0=initialtime;
    ocgTrj.timehorizon=initialendtime;
    ocgTrj.modelparameter=parametervalue(ocObj);
    ocgTrj.modelname=modelname(ocObj);
    ocgTrj.type=objoctype;
    ocgTrj.objectivevalue=OCGRADSOL.salvagevalue(initialendtime,initstate,ocgTrj.modelparameter);
    ocgTrj=ocgradtrajectory(ocgTrj);

end
cst_y=costate(ocgTrj);
y=state(ocgTrj);
OCGRADCONT.initialstate=y(:,1);
OCGRADCONT.endcostate=cst_y(:,end);

% number of states and controls
[control_num,state_num]=numberofvariables(ocObj);
fn=fieldnames(control_num);
laststatecoord=0;
for ii=1:length(fn)
    if ~isempty(state_num.(fn{ii}))
        OCGRADCONT.statecoordinate.(fn{ii})=1:state_num.(fn{ii});
        laststatecoord=laststatecoord+state_num.(fn{ii});
    end
end
OCGRADCONT.objectivecoordinate=laststatecoord+1;
OCGRADCONT.control_num=control_num;
OCGRADCONT.state_num=state_num;

%OCGRADCONT.statecoordinate=1:state_num.concentrated;
%OCGRADCONT.costatecoordinate=state_num.concentrated+(1:state_num.concentrated);
%OCGRADCONT.objectivecoordinate=OCGRADCONT.costatecoordinate(end)+1;
if isempty(fixinitialstate) % default: initial states are fixed and end states are free
    fixinitialstate=1:state_num.concentrated;
end
freeinitialstate=setdiff(1:state_num.concentrated,fixinitialstate);
freeendstate=setdiff(1:state_num.concentrated,fixendstate);

% information on control, mixed-control constraint
cstrnum=controlconstraintnum(ocObj);
mixedcstr=mixedconstraint(ocObj);

OCGRADCONT.cconstraint=0;
OCGRADCONT.mcconstraint=0;
if cstrnum
    OCGRADCONT.cconstraint=1;
    if mixedcstr
        OCGRADCONT.mcconstraint=1;
    else
        OCGRADCONT.mcconstraint=0;
    end
end
% filenames for storing intermediate continuation results and possible
% animation file
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

%
OCGRADCONT.basicmoviefilename=fullocmatfile(pathname,'AnimCont');

OCGRADCONT.modelname=modelname(ocObj);
OCGRADCONT.modelparameter=parametervalue(ocObj);

if ~isempty(modelname(ocgTrj)) && ~strcmp(modelname(ocObj),modelname(ocgTrj))
    ocmatmsg('Warning, model name of extremal and model are different')
end

% determine the free parameters
OCGRADCONT.method=method;
freeparameterctr=0;
sol.parameters=[];
% continuation type:  time,
OCGRADCONT.conttype=conttype;
switch conttype
    case 'time'
        if isempty(contindex)
            OCGRADCONT.contindex=2; % continue endtime
        else
            OCGRADCONT.contindex=contindex; % 1 initial time, 2 endtime
        end
    case 'initialstate'
        if isempty(targetvalue)
            return
        end
        if length(targetvalue)~=length(contindex)
            return
        end
    case 'endstate'
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
        OCGRADCONT.contindex=contindex;
    otherwise
        return
end
OCGRADCONT.targetvalue=targetvalue(:);
OCGRADCONT.octype=octype(ocgTrj);

sol.extremal=generategradextremal(ocgTrj);
switch conttype
    case 'time'
        if OCGRADCONT.contindex==1
            %sol.parameters=[sol.parameters initialtime(ocgTrj)];
            continuationparametervalue=initialtime(ocgTrj)+initcontpar;
        elseif OCGRADCONT.contindex==2
            %sol.parameters=[sol.parameters timehorizon(ocgTrj)+initcontpar];
            sol.extremal.t=(sol.extremal.t-sol.extremal.t0)/sol.extremal.t(end);
            continuationparametervalue=timehorizon(ocgTrj)+initcontpar;
        else
            return
        end
    case 'initialstate'
        %sol.parameters=0;
        continuationparametervalue=initcontpar;
        fixinitialstate=contindex;
        OCGRADCONT.initialstatecoordinate=contindex;
        OCGRADCONT.targetvector=OCGRADCONT.targetvalue-OCGRADCONT.initialstate(contindex);
    case 'endstate'
        %sol.parameters=0;
        continuationparametervalue=initcontpar;
        fixendstate=contindex;
        OCGRADCONT.endstate=y(:,end);

        OCGRADCONT.endstatecoordinate=contindex;
        OCGRADCONT.targetvector=OCGRADCONT.targetvalue-OCGRADCONT.endstate(contindex);
    case 'parameter'
        %sol.parameters=[sol.parameters parametervalue(ocObj,contindex)];
        continuationparametervalue=parametervalue(ocObj,contindex)+initcontpar;
        if ~isempty(OCGRADCONT.targetvalue)
            %    OCGRADCONT.targetvector=OCGRADCONT.targetvalue-OCGRADCONT.startingparametervalue;
        else
            %    OCGRADCONT.targetvector=1;
        end
end

% coordinate is used as a signifier for the position (row) in the solution
% path
% index is uses as a signifier for the position in the free parameters
% vector
OCGRADCONT.fixinitialstatecoordinate=fixinitialstate;
OCGRADCONT.freeinitialstatecoordinate=freeinitialstate;
OCGRADCONT.fixendstatecoordinate=fixendstate;
OCGRADCONT.freeendstatecoordinate=freeendstate;

OCGRADCONT.fixinitialstate=y(OCGRADCONT.fixinitialstatecoordinate,1);
OCGRADCONT.fixendstate=y(OCGRADCONT.fixendstatecoordinate,end);

if ~isempty(fixendstate)
    OCGRADCONT.freeendcostatecoordinate=fixendstate;
    freeparameterctr=freeparameterctr+(1:length(fixendstate));
    sol.parameters=[sol.parameters cst_y(fixendstate,end).'];
    OCGRADCONT.freeendcostateindex=fixendstate;
end

if ~isempty(freeinitialstate)
    OCGRADCONT.fixinitialcostatecoordinate=freeinitialstate;
    freeparameterctr=freeparameterctr+(1:length(fixendstate));
    sol.parameters=[sol.parameters cst_y(freeinitialstate,1).'];
    OCGRADCONT.fixinitialcostateindex=freeinitialstate;
end
OCGRADCONT.continuationindex=freeparameterctr(end)+1;
sol.parameters=[sol.parameters continuationparametervalue];

switch OCGRADCONT.octype
    case 'concentrated'
        n=length(time(ocgTrj));
        OCGRADCONT.loccontrolcoordinate=reshape(1:control_num.concentrated*n,[],n);
        OCGRADCONT.initiallocstatecoordinate=OCGRADCONT.loccontrolcoordinate(end)+(1:state_num.concentrated).';
        OCGRADCONT.endloccostatecoordinate=OCGRADCONT.initiallocstatecoordinate(end)+(1:state_num.concentrated).';
end
OCGRADCONT.TIMEMESH.num=length(sol.extremal.t);

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

