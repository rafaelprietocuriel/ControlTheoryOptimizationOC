function [xlc0,vlc0]= ocmatinit_H_LC(ocObj,dynPrim,ap,varargin)
%
% INIT_H_LC initialization of limit cycle continuation from Hopf
% bifurccation
%
% [X0,V0]= INIT_H_LC(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global lds OCMATFINITCONT 


xlc0=[];
vlc0=[];
if isempty(ocObj) || isempty(dynPrim)
    return
end

if ~isdynprimitive(dynPrim)
    error('Input argument is not a dynprimitve.')
end

if ~isequilibrium(dynPrim)
    error('Input argument is not an equilibrium.')
end
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end

OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
switch OCMATFINITCONT.modeltype
    case 'standardmodel'
        OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'CanonicalSystem');
        %OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'CanonicalSystemJacobian');
        %OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'CanonicalSystemParameterJacobian');
        OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
        OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
        OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
        OCMATFINITCONT.arcarg=arcargument(dynPrim);
    case 'odemodel'
        OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'Dynamics');
        OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'DynamicsJacobian');
        OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'DynamicsParameterJacobian');
        OCMATFINITCONT.arcarg=arcargument(dynPrim);
end
OCMATFINITCONT.arcarg=arcargument(dynPrim);

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');

[xlc0,vlc0]=init_H_LC(@modellimitcycle,dependentvar(dynPrim),parametervalue(ocObj).',ap,varargin{:});
