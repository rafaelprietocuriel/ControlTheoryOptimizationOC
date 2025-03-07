function [x0,v0]=ocmatinit_BP_EP(ocObj,dynPrim,s,h,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_H_H(OCOBJ,dynPrim,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT eds


x0=[];
v0=[];
idx=[];
userfunc='';
plotcoord=[];
if isempty(ocObj) || isempty(dynPrim)
    return
end
if isa(dynPrim,'gdynprimitive')
    [x0,v0]=ocmatinit_BP_EP_implicit(dynPrim,s,h,varargin);
    return
end

if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
if nargin>=5
    userfunc=varargin{1};
end
if nargin>=6
    idx=varargin{2};
end
if nargin>=7
    plotcoord=varargin{3};
end
if nargin>=8
    v0=varargin{4};
end
depvar=dependentvar(dynPrim);
if isempty(idx)
    idx=1:length(depvar);
end
if isempty(plotcoord)
    plotcoord=1;
end
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
OCMATFINITCONT.parametername=parametername(ocObj,eds.ActiveParams);
OCMATFINITCONT.conttype='equilibrium';

if isempty(userfunc)
    switch OCMATFINITCONT.modeltype
        case {'standardmodel','differentialgame'}
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquation');
            OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
            OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
            OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
            OCMATFINITCONT.plotcont=modelspecificfunc(ocObj,'PlotEquilibriumContinuation');
            
            OCMATFINITCONT.arcarg=arcargument(dynPrim);
        case 'odemodel'
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'Dynamics');
            OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'DynamicsJacobian');
            OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'DynamicsParameterJacobian');
            OCMATFINITCONT.plotcont=modelspecificfunc(ocObj,'PlotEquilibriumContinuation');

            OCMATFINITCONT.arcarg=arcargument(dynPrim);
    end
else
    funch=feval(modelspecificfunc(ocObj,userfunc));
    OCMATFINITCONT.dynamics=funch{1};
    OCMATFINITCONT.jacobian=funch{2};
    OCMATFINITCONT.parameterjacobian=funch{3};
    OCMATFINITCONT.arcarg=arcargument(dynPrim);
end

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');
OCMATFINITCONT.plotcoord=plotcoord;
OCMATFINITCONT.activeparameter=eds.ActiveParams;

try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case {'standardmodel','differentialgame'}
            [x0,v0]=init_BP_EP(@modelequilibrium,dependentvar(dynPrim),parametervalue(ocObj).',s,h,varargin{:});
            %[x0,v0]=init_BP_BP(@modelequilibrium,dependentvar(dynPrim),parametervalue(ocObj).',s,h,varargin{:});
            %[x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
        case 'odemodel'
            [x0,v0]=init_BP_EP(@odemodelequilibrium,dependentvar(dynPrim),parametervalue(ocObj).',s,h,varargin{:});
    end
    
catch
    lasterr
    ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end