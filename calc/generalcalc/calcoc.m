function sout=calcoc(varargin)

global OCMATCALC 
[OCMATCALC.problem_func,u0,v0,w0,opt]=ParseCommandLine(varargin{:});

if ~isfield(OCMATCALC,'xmesh')
    OCMATCALC.xmesh=[];
end

% function handles of the actual continuation type
problemhandles=feval(OCMATCALC.problem_func);
switch varargin{1}
    case 'godesol'
        OCMATCALC.funS=problemhandles{1};
    case 'gnlodesol'
        OCMATCALC.funF=problemhandles{1};
end
OCMATCALC.fun_adj=problemhandles{2};
OCMATCALC.funGradHam=problemhandles{3};
OCMATCALC.funJ=problemhandles{4};

if OCMATCALC.Nq>0
    OCMATCALC.funPQ=problemhandles{7};
    OCMATCALC.funEtaZeta=problemhandles{8};
end
if OCMATCALC.Ny>0
    if ~isfield(OCMATCALC,'xiT') ||isempty(OCMATCALC.xiT)
        OCMATCALC.funl_y=problemhandles{8};
    else
        OCMATCALC.funl_y=OCMATCALC.xiT;
    end
end
if OCMATCALC.Nz>0
    if ~isfield(OCMATCALC,'thetaT') ||isempty(OCMATCALC.thetaT)
        OCMATCALC.funl_z=problemhandles{6};
    else
        OCMATCALC.funl_z=OCMATCALC.thetaT;
    end
end
if isfield(OCMATCALC,'y0') && OCMATCALC.Ny>0
    if isempty(OCMATCALC.y0)
        OCMATCALC.funY0=problemhandles{9};
    else
        OCMATCALC.funY0=OCMATCALC.y0;
    end
end
if isfield(OCMATCALC,'z0') && OCMATCALC.Nz>0
    if isempty(OCMATCALC.z0)
        OCMATCALC.funZ0=problemhandles{5};
    else
        OCMATCALC.funZ0=OCMATCALC.z0;
    end
end
OCMATCALC.plotcalcoc=problemhandles{11};
OCMATCALC.printcontinuation=problemhandles{12};
OCMATCALC.saveintermediate=problemhandles{13};
OCMATCALC.formatsolution=problemhandles{22};

%opt=OCMATCALC.options(opt);
switch OCMATCALC.optimizationtype
    case 'min'
        lipschitzsign=1;
    case 'max'
        lipschitzsign=-1;
    otherwise
end
OCMATCALC.options.ConditionNumber=getocoptions(opt,'GRADIENT','ConditionNumber');
OCMATCALC.options.LipschitzConst=lipschitzsign*abs(getocoptions(opt,'GRADIENT','LipschitzConst'));
OCMATCALC.options.MaxIter=getocoptions(opt,'GRADIENT','MaxIter');
OCMATCALC.options.Display=getocoptions(opt,'GRADIENT','Display');
OCMATCALC.options.LineSrch=getocoptions(opt,'GRADIENT','LineSrch');

if ~isempty(getocoptions(opt,'GRADIENT','Projection'))
    OCMATCALC.options.Projection=getocoptions(opt,'GRADIENT','Projection');
end

OCMATCALC.options.UpperBoundV=getocoptions(opt,'GRADIENT','UpperBoundV');
OCMATCALC.options.LowerBoundV=getocoptions(opt,'GRADIENT','LowerBoundV');
OCMATCALC.options.UpperBoundU=getocoptions(opt,'GRADIENT','UpperBoundU');
OCMATCALC.options.LowerBoundU=getocoptions(opt,'GRADIENT','LowerBoundU');
OCMATCALC.GradSolver=getocoptions(opt,'GRADIENT','Solver');
if ischar(OCMATCALC.GradSolver)
    OCMATCALC.GradSolver=str2func(OCMATCALC.GradSolver);
end
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');
PlotCalc=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');

OCMATCALC.SaveIntermediate=SaveIntermediate;

if PlotCalc
    clf
end
StartTime = clock;

% initialize user workspace
if WorkSpace
    if OCMATCALC.probleminit(tmesh0,coeff0,tangent0)~=0
        ocmaterror('Initializer failed.');
    end
end

s.index=1;
s.label='00';
%s.data.sol=OCMATCALC.formatsolution(tmesh0,coeff0,tangent0);
s.msg='This is the first solution of the BVP continuation';


%sout=s;
[u,v,w,Val,du,dv,dw,y,z,p,q,xi,theta,eta,zeta]=OCMATCALC.GradSolver(u0,v0,w0,OCMATCALC.gradmap,OCMATCALC.options,OCMATCALC);
sout=OCMATCALC.formatsolution(OCMATCALC.xmesh,OCMATCALC.tspan,y,z,p,q,xi,theta,eta,zeta,u,v,w,Val,du,dv,dw);

if SaveIntermediate
    failed=OCMATCALC.saveintermediate(sout);
end
if PlotCalc
    PlotResult(u,v,w,Val,du,dv,dw,y,z,p,q,xi,theta,eta,zeta);
end


fprintf('\n');
EndTime = clock;
fprintf('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));

%---------------------------------
%
% Plot the solution
%
%
function b=PlotResult(u,v,w,Val,du,dv,dw,y,z,p,q,xi,theta,eta,zeta)
global OCMATCALC
b=OCMATCALC.plotcalcoc(OCMATCALC.xmesh,OCMATCALC.tspan,y,z,p,q,xi,theta,eta,zeta,u,v,w,Val,du,dv,dw);


%-----------------------------------
%
% Command line parser
%
%-----------------------------------

function [problemfunc,u,v,w,opt]=ParseCommandLine(problemfunc,solinit,varargin)

if nargin < 2
    ocmaterror('wrong number of input arguments');
end
opt=[];
if nargin > 2
    opt=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isfield(solinit,'u')
    u=solinit.u;
else
    u=[];
end
if isfield(solinit,'v')
    v=solinit.v;
else
    v=[];
end
if isfield(solinit,'w')
    w=solinit.w;
else
    w=[];
end
%--< END OF CMDL PARSER >--


%----------------------------
%
% DefaultProcessor
%
%----------------------------

function [failed,f,s]=DefaultProcessor(tmesh,coeff,tangent,s)
global OCMATCALC
% WM: this now actually calls the default processor,
% either with or without a singular point structure

if nargin > 3
    [failed,f,s]=OCMATCALC.defaultprocessor(tmesh,coeff,tangent,s);
else
    [failed,f]=OCMATCALC.defaultprocessor(tmesh,coeff,tangent);
end

