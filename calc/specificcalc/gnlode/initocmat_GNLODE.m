function sol=initocmat_GNLODE(ocObj,varargin)
%
% INITOCMAT_GNLODE initialization for general ODE problem

clear global OCMATCALC OCMATGNLODE
global OCMATCALC OCMATGNLODE

if nargin==1
    sol=[];
    return
end
y0=[];
nt=[];
nx=[];
t=[];
x=[];
z0=[];
q=[];
u=[];
v=[];
xiT=[];
thetaT=[];
startidx=1;
if isstruct(varargin{1})
    t=varargin{1}.t;
    x=varargin{1}.x;
    u=varargin{1}.u;
    v=varargin{1}.v;
    w=varargin{1}.w;
    try
        y0=varargin{1}.modelinfo.y0;
    end
    try
        z0=varargin{1}.modelinfo.z0;
    end
    nt=length(t);
    nx=length(x);
    startidx=2;
end
if mod(nargin-startidx,2)
    error('Input arguments have to come in pairs or as structure.')
end
for ii=startidx:2:nargin-1
    if ~ischar(varargin{ii})
        error('Characterization has to be a string variable.')
    end
    switch lower(varargin{ii})
        case {'linitialstate','y0'}
            y0=varargin{ii+1};
        case {'nlinitialstate','z0'}
            z0=varargin{ii+1};
        case {'linitialstate','xit'}
            xiT=varargin{ii+1};
        case {'nlinitialstate','thetat'}
            thetaT=varargin{ii+1};
        case {'nlcontrol','u'}
            u=varargin{ii+1};
        case {'lcontrol','v'}
            v=varargin{ii+1};
        case {'lintvar','q'}
            q=varargin{ii+1};
        case {'time','t'}
            t=varargin{ii+1};
        case {'parameter','space','x'}
            x=varargin{ii+1};
        case 'nt'
            nt=varargin{ii+1};
        case 'nx'
            nx=varargin{ii+1};
        otherwise
    end
end
if length(t)==2
    t=linspace(t(1),t(2),nt);
end
if length(t)<3
    error('Time discretization consist of less than 3 points.')
end
if length(x)==2
    x=linspace(x(1),x(2),nx);
end
if length(x)<3
    error('Space/parameter discretization consist of less than 3 points.')
end

% initialize global variable (OCMATCONT) for general continuation process
OCMATCALC.modelname=modelname(ocObj);
OCMATCALC.modelfunc=modelspecificfunc(ocObj,'4GNLOdeCalc');

funch=OCMATCALC.modelfunc();

OCMATGNLODE.statedynamics=funch{1};
OCMATGNLODE.adjointdynamics=funch{2};
OCMATGNLODE.gradienthamiltonian=funch{3};
OCMATGNLODE.objectivefunction=funch{4};
OCMATGNLODE.initialstate=funch{5};
OCMATGNLODE.transversalitycondition=funch{6};
OCMATGNLODE.integralconstraint=funch{7};
OCMATGNLODE.adjointintegral=funch{8};

% general function
OCMATGNLODE.plotcalcoc=funch{11};
OCMATGNLODE.datapath=funch{12};
OCMATGNLODE.saveintermediatefiles=funch{13};


pathname=OCMATGNLODE.datapath();
[resultfile,globalvarfile]=OCMATGNLODE.saveintermediatefiles();
OCMATGNLODE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATGNLODE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATGNLODE.dt=diff(t);
OCMATGNLODE.dx=diff(x);
OCMATGNLODE.initialstate='fix';
OCMATGNLODE.z0=z0;

OCMATGNLODE.modelparameter=parametervalue(ocObj);
OCMATGNLODE.exp_rt=exp(-discountrate(ocObj)*t(:));

OCMATCALC.tspan = t.';
OCMATCALC.xmesh =x(:);
OCMATCALC.dt = diff(t);
OCMATCALC.Nt = nt;
OCMATCALC.Nx = nx;
OCMATCALC.Ny = nonlocalodestatenum(ocObj);
OCMATCALC.Nz = localodestatenum(ocObj);
OCMATCALC.Nu = nonlocalcontrolnum(ocObj);
OCMATCALC.Np = nonlocalintegralstatenum(ocObj);
OCMATCALC.Nq = localintegralstatenum(ocObj);
OCMATCALC.Nv = localcontrolnum(ocObj);
OCMATCALC.Nw = 0;
OCMATCALC.gradmap=@fun_gradmap_gnlodeoc;
OCMATCALC.optimizationtype=optimizationtype(ocObj);
OCMATCALC.modeltype=modeltype(ocObj);

OCMATCALC.z0=z0;
OCMATCALC.y0=y0;
OCMATCALC.xiT=xiT;
OCMATCALC.thetaT=thetaT;

sol.t=t;
sol.x=x;
sol.v=v;
sol.u=u;


% check consistency
if ~isempty(y0)
    [row cols]=size(y0);
    if row~=OCMATCALC.Nx || cols~=OCMATCALC.Ny
        error('Space/parameter distribution of initial non-local states are not consistent with discretization.')
    end
end
if ~isempty(z0)
    [row cols]=size(z0);
    if row~=OCMATCALC.Nz || cols~=OCMATCALC.Nz
        error('Space/parameter distribution of initial local states are not consistent with discretization.')
    end
end
if ~isempty(u)
    [row cols planes]=size(u);
    if row~=OCMATCALC.Nx || cols~=OCMATCALC.Nu || planes~=nt
        error('Space/parameter distribution of initial non-local controls are not consistent with discretization.')
    end
end
if ~isempty(v)
    [row cols]=size(v);
    if row~=OCMATCALC.Nv || cols~=OCMATCALC.Nt
        error('Space/parameter distribution of initial local controls are not consistent with discretization.')
    end
end
