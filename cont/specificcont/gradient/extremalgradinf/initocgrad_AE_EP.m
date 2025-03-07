function sol=initocgrad_AE_EP(ocObj,ocgAsym,conttype,contindex,varargin)
clear global OCGRADSOL OCGRADCONT
global OCGRADSOL OCGRADCONT
sol=[];
targetvalue=[];
targettype='';
opt=[];
truncationtime=[];
objoctype='';
initcontpar=[];
dist2ep=[];
testdistance=[];

targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targettypeidx=find(strcmpi(varargin,'targettype')); % contpar (default case), norm2eq (distance of last solutionpoint to the equilibrium)
truncationtimeidx=find(strcmpi(varargin,'truncationtime'));
octypeidx=find(strcmpi(varargin,'octype'));
optionidx=find(strcmpi(varargin,'option'));
initcontparidx=find(strcmpi(varargin,'initcontpar'));
dist2epidx=find(strcmpi(varargin,'dist2ep'));
testdistanceidx=find(strcmpi(varargin,'testdistance'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(octypeidx)
    objoctype=varargin{octypeidx+1}; % only used if ocgTrj is empty
end
if ~isempty(truncationtimeidx)
    truncationtime=varargin{truncationtimeidx+1};
end
if ~isempty(octypeidx)
    objoctype=varargin{octypeidx+1}; % only used if ocgEP is empty
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1}; % only used if ocgEP is empty
end
if ~isempty(initcontparidx)
    initcontpar=varargin{initcontparidx+1}; % only used if ocgEP is empty
end
if ~isempty(dist2epidx)
    dist2ep=varargin{dist2epidx+1};
end
if ~isempty(testdistanceidx)
    testdistance=varargin{testdistanceidx+1};
end
if isempty(objoctype)
    objoctype='concentrated'; % default is local (not distributed)
end
if isempty(initcontpar)
    initcontpar=0;
end
if isempty(targettype)
    targettype='contpar';
end
if isempty(dist2ep)
    dist2ep=0;
end
if isempty(testdistance)
    testdistance=0;
end
if isempty(opt)
    opt=defaultocoptions;
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

if dist2ep || strcmp(conttype,'time')
    OCGRADCONT.normalizetime=1;
else
    OCGRADCONT.normalizetime=0;
end
if isocgradasymptotic(ocgAsym)
    ocgEP=limitset(ocgAsym);
    ocgTrj=ocgradtrajectory(ocgAsym);
    truncationtime=timehorizon(ocgTrj);
    if abs(truncationtime-1)<1e-10
        ocmatmsg('Warning, check if truncation time 1 is correct.')
    end
    cst_y=costate(ocgTrj);
    cst_y=cst_y(:,end);
elseif isequilibrium(ocgAsym)
    if isempty(truncationtime)
        return
    end

    ocgEP=ocgAsym;
    cst_y=costate(ocgEP);
    cst_y=cst_y(:);
    % generate initial path
    if isempty(ocgEP)
        return
    end
    initocPt.y=dependentvar(ocgEP);
    initocPt.x=0;
    initocPt.arcarg=arcargument(ocgEP);
    initocPt.arcinterval=[0 0];
    initocPt.arcposition=[1;1];
    ctrl=control(ocgEP);
    %ctrl=control(ocObj,octrajectory(initocPt));

    % create initial ocgradtrajectory
    trivialarcmeshnum=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
    ocgTrj.argument.t=linspace(0,truncationtime,trivialarcmeshnum);
    ocgTrj.variable.y=state(ocgEP);
    ocgTrj.variable.y=ocgTrj.variable.y(:,ones(1,trivialarcmeshnum));
    ocgTrj.variable.cst_y=costate(ocgEP);
    ocgTrj.variable.cst_y=ocgTrj.variable.cst_y(:,ones(1,trivialarcmeshnum));
    ocgTrj.variable.v=ctrl;
    ocgTrj.variable.v=ocgTrj.variable.v(:,ones(1,trivialarcmeshnum));
    ocgTrj.variable.o=0;
    ocgTrj.variable.o=ocgTrj.variable.o(:,ones(1,trivialarcmeshnum));
    ocgTrj.initargument.t0=0;
    ocgTrj.timehorizon=truncationtime;
    ocgTrj.modelparameter=parametervalue(ocObj);
    ocgTrj.modelname=modelname(ocObj);
    ocgTrj.type=objoctype;
    ocgTrj.objectivevalue=0;
    ocgTrj=ocgradtrajectory(ocgTrj);
else
    return
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
OCGRADSOL.asymptotictransversaltycondition=funch{15}{1}; % asymptotic transversality condition
OCGRADSOL.asymptoticobjectivevalue=funch{15}{2}; % objective value H/r
OCGRADSOL.asymptoticobjectivefunction=funch{15}{3}; % objective value H/r

OCGRADSOL.plotcontinuation=funch{11};
OCGRADSOL.datapath=funch{20};
OCGRADSOL.saveintermediatefiles=funch{21};

OCGRADSOL.equilibrium=ocgEP;

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

%%% generate initial octrajectory if not provided

% filenames for storing intermediate continuation results
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCGRADCONT.modelname=modelname(ocObj);
OCGRADCONT.modelparameter=parametervalue(ocObj);

OCGRADSOL.saddlepoint=dependentvar(ocgEP);
OCGRADSOL.Jacobian=linearization(ocgEP);
%OCGRADSOL.arcarg4equilibrium=arcargument(ocgEP);

[OCGRADSOL.asymptoticmatrix]=asymptoticbc(OCGRADSOL.Jacobian,'s','c',ZeroDeviationTolerance);
%
OCGRADCONT.initialstatevalue=state(ocgTrj);

% continuation type:  state, (end)time,
OCGRADCONT.conttype=conttype;
switch conttype
    case 'state'
        OCGRADCONT.initialstatevalue=OCGRADCONT.initialstatevalue(contindex,1);
        if isempty(targetvalue)
            ocmatmsg('Target value is empty.')
            return
        end
        OCGRADCONT.contvector=targetvalue(:)-OCGRADCONT.initialstatevalue;
        contparameter=initcontpar;
        sol.extremal=generategradextremal(ocgTrj,OCGRADCONT.normalizetime);
    case 'time'
        if isempty(contindex)
            contindex=2;
        end
        if contindex==2
            sol.extremal=generategradextremal(ocgTrj,OCGRADCONT.normalizetime); % normalize time to interval [0,1]
            contparameter=truncationtime;
        else
            return
        end
    otherwise
        return
end
sol.parameters=cst_y.';
OCGRADCONT.freeparametercoordinate=1:length(sol.parameters);
OCGRADCONT.dist2ep=dist2ep;
if dist2ep
    y=state(ocgTrj);
    cst_y=costate(ocgTrj);
    OCGRADCONT.dist2epvalue=sqrt(sum((OCGRADSOL.saddlepoint-[y(:,end);cst_y(:,end)]).^2));
    OCGRADCONT.normalizetime=1;
    sol.parameters=[sol.parameters truncationtime];
    OCGRADCONT.freetimehorizoncoordinate=length(sol.parameters);
end
sol.parameters=[sol.parameters contparameter];
OCGRADCONT.continuationcoordinate=length(sol.parameters);
if strcmp(conttype,'time') && contindex==2
    OCGRADCONT.freetimehorizoncoordinate=OCGRADCONT.continuationcoordinate;
end
[control_num,state_num]=numberofvariables(ocObj);

fn=fieldnames(control_num);
laststatecoord=0;
for ii=1:length(fn)
    if ~isempty(state_num.(fn{ii}))
        OCGRADCONT.statecoordinate.(fn{ii})=1:state_num.(fn{ii});
        OCGRADCONT.costatecoordinate.(fn{ii})=OCGRADCONT.statecoordinate.(fn{ii})(end)+(1:state_num.(fn{ii}));
        laststatecoord=laststatecoord+state_num.(fn{ii});
    end
end
OCGRADCONT.objectivecoordinate=laststatecoord+1;
OCGRADCONT.control_num=control_num;
OCGRADCONT.state_num=state_num;
OCGRADCONT.targetvalue=targetvalue(:);
OCGRADCONT.targettype=targettype;
OCGRADCONT.octype=octype(ocgTrj);
OCGRADCONT.contindex=contindex;
OCGRADCONT.testdistance=testdistance;

switch OCGRADCONT.octype
    case 'concentrated'
        n=length(time(ocgTrj));
        OCGRADCONT.loccontrolcoordinate=reshape(1:control_num.concentrated*n,[],n);
        OCGRADCONT.initiallocstatecoordinate=OCGRADCONT.loccontrolcoordinate(end)+(1:state_num.concentrated).';
        OCGRADCONT.endloccostatecoordinate=OCGRADCONT.initiallocstatecoordinate(end)+(1:state_num.concentrated).';
end
OCGRADCONT.TIMEMESH.num=length(sol.extremal.t);

%
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