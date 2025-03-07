function sol=initocmat_ODE_CONT(ocObj,solinit,targetparametervalue,contfile,varargin)
% INITOCMAT_ODE_CONT initialization for
%
% SOL=INITOCMAT_ODE_CONT(OCOBJ,SOLINIT,TARGETVALUE)


clear global OCMATCONT OCMATODEBVP
global OCMATCONT OCMATODEBVP

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
OCMATODEBVP.plotcontinuation=funch{11};
OCMATODEBVP.testadmissibility=funch{12};
OCMATODEBVP.datapath=funch{20};
OCMATODEBVP.saveintermediatefiles=funch{21};
OCMATODEBVP.generatesolstruct=funch{22};

sol=generatesolstruct(solinit);
sol=OCMATODEBVP.generatesolstruct(sol);
% mode and path specific variables
OCMATODEBVP.parametervalue=parametervalue(ocObj);
OCMATODEBVP.initialtime=sol.x(1);
OCMATODEBVP.arcinterval=sol.arcinterval;
OCMATODEBVP.targetparametervalue=targetparametervalue;
OCMATCONT.codimension=1;

pathname=OCMATODEBVP.datapath();
[resultfile,globalvarfile]=OCMATODEBVP.saveintermediatefiles();
OCMATODEBVP.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATODEBVP.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

function sol=generatesolstruct(ocTrj)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];

if isfield(ocTrj.solverinfo,'tangent')
    sol.solverinfo.tangent=ocTrj.solverinfo.tangent;
else
    sol.solverinfo.tangent=[];
end
if isfield(ocTrj.solverinfo,'coeff')
    sol.solverinfo.coeff=ocTrj.solverinfo.coeff;
else
    sol.solverinfo.coeff=[];
end
if isfield(ocTrj.solverinfo,'yp')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ypmid')
    sol.solverinfo.ypmid=ocTrj.solverinfo.ypmid;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.ymid=ocTrj.solverinfo.ymid;
end