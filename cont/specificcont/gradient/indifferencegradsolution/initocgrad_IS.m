function sol=initocgrad_IS(ocObj,ocgMTrj,conttype,contindex,varargin)
clear global OCGRADSOL OCGRADCONT
global OCGRADSOL OCGRADCONT
sol=[];
freeparameter=[];
targetvalue=[];
targetindex=[];
freestatevector=[];
freeparameterindex=[];
opt=[];

if ischar(contindex)
    contindex=parameterindex(ocObj,contindex);
end
freeparameteridx=find(strcmpi(varargin,'freeparameter')); %
freestatevectoridx=find(strcmpi(varargin,'freestatevector')); %
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetindexidx=find(strcmpi(varargin,'targetindex'));
optionidx=find(strcmpi(varargin,'option'));

if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
    if ischar(freeparameter)
        freeparameterindex=parameterindex(ocObj,freeparameter);
    else
        freeparameterindex=freeparameter;
        freeparameter=parametername(ocObj,freeparameterindex);
        freeparameter=strtrim([char(freeparameter) char(repmat({','},length(freeparameterindex),1))]);
        freeparameter=strrep(reshape(freeparameter.',1,[]),' ','');
        freeparameter(end)=[];
        freeparameter=['[', freeparameter ']'];
    end
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targetindexidx)
    targetindex=varargin{targetindexidx+1};
end
if ~isempty(freestatevectoridx)
    freestatevector=varargin{freestatevectoridx+1};
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1}; % only used if ocEP is empty
end
if isempty(opt)
    opt=defaultocoptions;
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

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
if isocgradasymptotic(ocgMTrj)
    OCGRADSOL.asymptotictransversaltycondition=funch{15}{1};
    OCGRADSOL.asymptoticobjectivevalue=funch{15}{2}; % objective value H/r
    OCGRADSOL.asymptoticobjectivefunction=funch{15}{3}; % objective value H/r
end
OCGRADSOL.plotcontinuation=funch{11};
OCGRADSOL.datapath=funch{20};
OCGRADSOL.saveintermediatefiles=funch{21};

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

% filenames for storing intermediate continuation results
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCGRADCONT.modelname=modelname(ocObj);
OCGRADCONT.modelparameter=parametervalue(ocObj);

OCGRADSOL.degree=multiplicity(ocgMTrj);
OCGRADCONT.octype=octype(ocgMTrj);
OCGRADCONT.targetvalue=targetvalue;

[control_num,state_num]=numberofvariables(ocObj);
OCGRADCONT.continuationindex=OCGRADSOL.degree;
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
% continuation type:  time,
OCGRADCONT.conttype=conttype;
OCGRADCONT.contindex=contindex;
OCGRADCONT.numfreestatevector=[];

sol.parameters=[];
parametercounter=0;
for ii=1:OCGRADSOL.degree
    switch class(ocgMTrj(ii))
        case 'ocgradasymptotic'
            ocEP=limitset(ocgMTrj);
            OCGRADSOL.equilibrium{ii}=ocEP(ii);
            OCGRADSOL.saddlepoint{ii}=dependentvar(ocEP(ii));
            OCGRADSOL.Jacobian{ii}=linearization(ocEP(ii));
            [OCGRADSOL.asymptoticmatrix{ii}]=asymptoticbc(OCGRADSOL.Jacobian{ii},'s','c',ZeroDeviationTolerance);
            cst_y=costate(ocgMTrj(ii));
            sol.parameters=[sol.parameters cst_y(:,end).'];
            OCGRADCONT.freecostatecoordinate{ii}=(parametercounter+1):length(sol.parameters);
            OCGRADCONT.freecostatenum{ii}=length(OCGRADCONT.freecostatecoordinate{ii});
            parametercounter=length(sol.parameters);
            OCGRADCONT.trajectoryclass{ii}='inf';
        otherwise
            OCGRADCONT.trajectoryclass{ii}='fin';
    end

end
OCGRADCONT.numfreestatevector=0;
OCGRADCONT.freestatevector=freestatevector;
OCGRADCONT.freeparametercoordinate=[];
if size(freestatevector,2)==1
    OCGRADCONT.numfreestatevector=1;
elseif size(freestatevector,2)==2
    OCGRADCONT.numfreestatevector=2;
end
switch conttype
    case 'initialstate'
        contparametervalue=0;
    case 'parameter'
        contparametervalue=parametervalue(ocObj,contindex);
    otherwise
        return
end
if ~isempty(freeparameterindex)
    sol.parameters=[sol.parameters parametervalue(ocObj,freeparameterindex)];
    OCGRADCONT.freeparametercoordinate=(parametercounter+1):length(sol.parameters);
    parametercounter=length(sol.parameters);
elseif OCGRADCONT.numfreestatevector>=1
    sol.parameters=[sol.parameters 0];
    OCGRADCONT.freestatevectorcoordinate=(parametercounter+1):length(sol.parameters);
    parametercounter=length(sol.parameters);
end
OCGRADCONT.freeparameter=freeparameter;
OCGRADCONT.freeparameterindex=freeparameterindex;

OCGRADCONT.initialstatevalue=initialstate(ocgMTrj(1));

for ii=1:OCGRADSOL.degree
    sol.extremal(ii)=generategradextremal(ocgMTrj(ii));
end
sol.parameters=[sol.parameters contparametervalue];
OCGRADCONT.continuationcoordinate=length(sol.parameters);

if ischar(targetindex)
    targetindex=parameterindex(ocObj,targetindex);
end
OCGRADCONT.targetindex=targetindex;
if ~isempty(targetindex)
    if strcmp(conttype,'parameter')
        if targetindex==contindex
            OCGRADCONT.targetcoordinate=OCGRADCONT.continuationcoordinate;
        elseif any(targetindex==freeparameterindex)
            OCGRADCONT.targetcoordinate=OCGRADCONT.freeparametercoordinate(targetindex==freeparameterindex);
        end
    end
end
offset=0;
for ii=1:OCGRADSOL.degree
    switch OCGRADCONT.octype{ii}
        case 'concentrated'
            n=length(time(ocgMTrj(ii)));
            OCGRADCONT.loccontrolcoordinate{ii}=reshape(1:control_num.concentrated*n,[],n)+offset;
            OCGRADCONT.initiallocstatecoordinate{ii}=OCGRADCONT.loccontrolcoordinate{ii}(end)+(1:state_num.concentrated).';
            OCGRADCONT.endloccostatecoordinate{ii}=OCGRADCONT.initiallocstatecoordinate{ii}(end)+(1:state_num.concentrated).';
            offset=OCGRADCONT.endloccostatecoordinate{ii}(end);
    end
end
OCGRADCONT.TIMEMESH.num=length(time(ocgMTrj(1))); % assume that all indifferencesolutions have the same number of the time grid

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
