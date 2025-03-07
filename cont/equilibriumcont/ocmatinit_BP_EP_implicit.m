function [x0,v0]=ocmatinit_BP_EP_implicit(dynPrim,s,h,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_H_H(OCOBJ,dynPrim,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT eds


x0=[];
v0=[];
userfunc='';
plotcoord=[];
if isempty(dynPrim)
    return
end
ocObj=loadmodel(dynPrim);
if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
if nargin>=5
    userfunc=varargin{1};
end
if nargin>=7
    plotcoord=varargin{3};
end
if nargin>=8
    v0=varargin{4};
end
depvar=dependentvar(dynPrim);
depvar=depvar{1};
if isempty(plotcoord)
    plotcoord=1;
end
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
OCMATFINITCONT.parametername=parametername(ocObj,eds.ActiveParams);
OCMATFINITCONT.conttype='equilibrium';

if isempty(userfunc)
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquation');
            OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
            OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
            OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
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
        case 'standardmodel'
            [x0,v0]=init_BP_EP(@modelequilibrium,depvar,parametervalue(ocObj).',s,h);
            %[x0,v0]=init_BP_BP(@modelequilibrium,dependentvar(dynPrim),parametervalue(ocObj).',s,h,varargin{:});
            %[x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
    end
    
catch
    lasterr
    ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end