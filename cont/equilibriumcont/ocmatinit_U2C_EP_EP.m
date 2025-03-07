function [x0,v0]=ocmatinit_U2C_EP_EP(ocObj,dynPrim,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_EP_EP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];
cctrl='';
ccostate=[];
plotcoord=[];
arcarg=[];
if isempty(ocObj) || isempty(dynPrim)
    return
end
if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
ap=parameternum(ocObj)+1;
if nargin>=3
    cctrl=varargin{1};
end
if nargin>=4
    ccostate=varargin{2};
end
if nargin>=5
    plotcoord=varargin{3};
end
if nargin>=6
    v0=varargin{4};
end
if nargin>=7
    arcarg=varargin{5};
end
depvar=dependentvar(dynPrim);
n=statenum(ocObj);
m=controlnum(ocObj);

if isempty(cctrl)
    cctrl=zeros(m,1);
end
if isempty(ccostate)
    ccostate=zeros(n,1);
end
if isempty(plotcoord)
    plotcoord=2;
end
depvar=[depvar;ccostate];

OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
switch OCMATFINITCONT.modeltype
    case 'standardmodel'
        OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquationCont');
        OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
        OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
        OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
        OCMATFINITCONT.plotcont=modelspecificfunc(ocObj,'PlotEquilibriumContinuation');

        if isempty(arcarg)
            OCMATFINITCONT.arcarg=arcargument(dynPrim);
        else
            OCMATFINITCONT.arcarg=arcarg;
        end
        OCMATFINITCONT.cctrl=cctrl;
        OCMATFINITCONT.ccostate=ccostate;
    case 'odemodel'
        OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'Dynamics');
        OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'DynamicsJacobian');
        OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
        OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'DynamicsParameterJacobian');
        OCMATFINITCONT.plotcont=modelspecificfunc(ocObj,'PlotEquilibriumContinuation');
        OCMATFINITCONT.arcarg=arcargument(dynPrim);
end
modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');
OCMATFINITCONT.plotcoord=plotcoord;
try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            [x0,v0]=init_EP_EP(@modelequilibriumu2c,depvar,[parametervalue(ocObj) 0],ap);
            %[x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
        case 'odemodel'
            [x0,v0]=init_EP_EP(@modelequilibriumu2c,depvar,[parametervalue(ocObj) 0],ap,varargin{3:end});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end