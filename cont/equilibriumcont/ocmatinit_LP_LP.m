function [x0,v0]=ocmatinit_LP_LP(ocObj,dynPrim,ap,varargin)
%
% ocmatinit_LP_LP Initializes a fold continuation from a limitpoint.
%
% [X0,V0]= ocmatinit_LP_LP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];

if isempty(ocObj) || isempty(dynPrim)
    return
end
if isa(dynPrim,'gdynprimitive')
    [x0,v0]=ocmatinit_LP_LP_implicit(dynPrim,ap,varargin{:});
    return
end

if isnumeric(dynPrim)
    dynPrim.y=dynPrim;
    dynPrim.x=0;
    dynPrim.arcarg=0;
    dynPrim=dynprimitive(dynPrim);
elseif ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
if isempty(ap) || any(ap>parameternum(ocObj)) || any(ap<1)
    ocmaterror('Parameter index ''%s'' is not valid',char(ap))
end
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
switch OCMATFINITCONT.modeltype
    case 'standardmodel'
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
        OCMATFINITCONT.arcarg=arcargument(dynPrim);
end

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');

try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            [x0,v0]=init_LP_LP(@modellimitpoint,dependentvar(dynPrim),parametervalue(ocObj).',ap,varargin{:});
        case 'odemodel'
            [x0,v0]=init_LP_LP(@modellimitpoint,dependentvar(dynPrim),parametervalue(ocObj).',ap,varargin{:});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end