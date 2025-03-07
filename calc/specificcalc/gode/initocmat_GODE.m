function sol=initocmat_GODE(ocObj,varargin)
%
% INITOCMAT_GODE initialization for general ODE problem

clear global OCMATCALC OCMATGODE
global OCMATCALC OCMATGODE

if nargin==1
    sol=[];
    return
end
nt=[];
t=[];
z0=[];
q=[];
v=[];
thetaT=[];

nums=statenum(ocObj);
if isstruct(varargin{1})
    t=varargin{1}.t;
    z0=varargin{1}.z(1:nums,1);
    v=varargin{1}.v;
    nt=length(t);
else
    error(nargchk(4, inf, nargin))
    z0=varargin{1};
    v=varargin{2};
    if nargin==4
        t=varargin{3};
    else
        nt=varargin{4};
        if length(varargin{3})==1
            timeinterval=[0 varargin{3}];
        else
            timeinterval=varargin{3};
        end
        t=linspace(timeinterval(1),timeinterval(2),nt);
    end
end

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCALC.modelname=modelname(ocObj);
OCMATCALC.modelfunc=modelspecificfunc(ocObj,'4GOdeCalc');

funch=OCMATCALC.modelfunc();

OCMATGODE.statedynamics=funch{1};
OCMATGODE.adjointdynamics=funch{2};
OCMATGODE.gradienthamiltonian=funch{3};
OCMATGODE.objectivefunction=funch{4};
OCMATGODE.initialstate=funch{5};
OCMATGODE.transversalitycondition=funch{6};

% general function
OCMATGODE.plotcalcoc=funch{11};
OCMATGODE.datapath=funch{12};
OCMATGODE.saveintermediatefiles=funch{13};


pathname=OCMATGODE.datapath();
[resultfile,globalvarfile]=OCMATGODE.saveintermediatefiles();
OCMATGODE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATGODE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATGODE.dt=diff(t);
OCMATGODE.initialstate='fix';
OCMATGODE.z0=z0;

OCMATGODE.modelparameter=parametervalue(ocObj);
OCMATGODE.exp_rt=exp(-discountrate(ocObj)*t);

OCMATCALC.tspan = t.';
OCMATCALC.dt = diff(t);
OCMATCALC.Nt = nt;
OCMATCALC.Nz = nums;
OCMATCALC.Nv = controlnum(ocObj);
OCMATCALC.Np = 0;
OCMATCALC.Nq = 0;
OCMATCALC.Nu = 0;
OCMATCALC.Nw = 0;
OCMATCALC.gradmap=@fun_gradmap_gode;
OCMATCALC.optimizationtype=optimizationtype(ocObj);
sol.t=t;
sol.v=v;
