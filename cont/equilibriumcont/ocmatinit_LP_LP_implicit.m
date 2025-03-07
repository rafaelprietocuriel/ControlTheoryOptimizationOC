function [x0,v0]=ocmatinit_LP_LP_implicit(dynPrim,ap,varargin)
%
% ocmatinit_LP_LP Initializes a fold continuation from a limitpoint.
%
% [X0,V0]= ocmatinit_LP_LP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];

if isempty(dynPrim)
    return
end

ocObj=loadmodel(dynPrim);
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
if isempty(ap) || any(ap>parameternum(ocObj)) || any(ap<1)
    ocmaterror('Parameter index ''%s'' is not valid',char(ap))
end
depvar=dependentvar(dynPrim);
depvar=depvar{1};

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
end

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');

try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            [x0,v0]=init_LP_LP(@modellimitpoint,depvar,parametervalue(ocObj).',ap,varargin{:});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end