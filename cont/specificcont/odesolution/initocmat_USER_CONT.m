function sol=initocmat_USER_CONT(ocObj,solinit,contfile,parindex,varargin)
% INITOCMAT_ODE_CONT initialization for
%
% SOL=INITOCMAT_ODE_CONT(OCOBJ,SOLINIT,TARGETVALUE)


clear global OCMATCONT OCMATODEBVP
global OCMATCONT OCMATODEBVP

targetparametervalue=[];
endtime=[];
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
endtimeidx=find(strcmpi(varargin,'endtime'));

if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(endtimeidx)
    endtime=varargin{endtimeidx+1};
end
if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
% if isempty(parindex)
%     ocmatmsg('No or unknown parameter specified.')
%     return
% end
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,contfile);

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATODEBVP.dynamics=funch{1};
OCMATODEBVP.jacobian=funch{2}{1};
OCMATODEBVP.parameterjacobian=funch{2}{2};

% function for the boundary conditions
OCMATODEBVP.bc=funch{5};

% function for Jacobian
OCMATODEBVP.bcjacobian=funch{6};

% general function
OCMATODEBVP.probleminit=funch{10};
OCMATODEBVP.plotcontinuation=funch{11};
OCMATODEBVP.testadmissibility=funch{12};
OCMATODEBVP.datapath=funch{20};
OCMATODEBVP.saveintermediatefiles=funch{21};

OCMATCONT.DOMAINDDATA.numode=size(solinit.y,1);
OCMATCONT.DOMAINDDATA.eqcoord=1:size(solinit.y,1);
OCMATCONT.DOMAINDDATA.odecoord=1:size(solinit.y,1);
OCMATODEBVP.parametervalue=parametervalue(ocObj);
sol=generatesolstruct(solinit,endtime);
sol.parameters=OCMATODEBVP.parametervalue(parindex);

% mode and path specific variables
OCMATODEBVP.initialtime=sol.x(1);
OCMATODEBVP.switchtimecoord=[];
OCMATODEBVP.autonomous=1;
OCMATODEBVP.endtime=sol.arcinterval(end);
OCMATCONT.codimension=1;

OCMATODEBVP.varyparameterindex=parindex;
OCMATODEBVP.targetparametervalue=targetparametervalue;

pathname=OCMATODEBVP.datapath();
[resultfile,globalvarfile]=OCMATODEBVP.saveintermediatefiles();
OCMATODEBVP.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATODEBVP.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);


function sol=generatesolstruct(solinit,endtime)

if isoctrajectory(solinit)
    sol.x=independentvar(solinit);
    sol.y=dependentvar(solinit);
    sol.arcarg=arcargument(solinit);
    sol.x0=sol.x(1);
    sol.arcinterval=arcinterval(solinit);
    sol.arcposition=arcposition(solinit);
    sol.parameters=[];
    sol.solver='';%solvername;
else
    sol.x=solinit.x;
    sol.y=solinit.y;
    sol.arcarg=0;
    sol.x0=sol.x(1);
    sol.arcposition=[1;length(sol.x)];
    if ~isempty(endtime)
        sol.arcinterval=[0 endtime];
    end
    sol.parameters=[];
    sol.solver='';%solvername;
end