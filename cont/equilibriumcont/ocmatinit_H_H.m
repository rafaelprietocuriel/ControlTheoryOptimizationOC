function [x0,v0]=ocmatinit_H_H(ocObj,ocEP,ap,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_H_H(OCOBJ,OCEP,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];
if isempty(ocObj) || isempty(ocEP)
    return
end
OCMATFINITCONT.conttype='hopf';

if ~isequilibrium(ocEP)
    ocmaterror('Input argument is not an equilibrium.')
end
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
OCMATFINITCONT.parametername=parametername(ocObj,ap);
if isempty(ap) || max(ap)>parameternum(ocObj) || length(ap)<2
    ocmaterror('Parameter index ''%s'' is not valid',char(ap))
end
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquation');
OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
OCMATFINITCONT.testadmissibility=modelspecificfunc(ocObj,'EquilibriumAdmissible');
OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
OCMATFINITCONT.hessian=modelspecificfunc(ocObj,'CanonicalSystemHessian');
OCMATFINITCONT.parameterhessian=[];%modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
OCMATFINITCONT.der3=modelspecificfunc(ocObj,'CanonicalSystemDer3');
OCMATFINITCONT.der4=modelspecificfunc(ocObj,'CanonicalSystemDer4');
OCMATFINITCONT.der5=modelspecificfunc(ocObj,'CanonicalSystemDer5');
OCMATFINITCONT.arcarg=arcargument(ocEP);

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');

try
    % calls MATCONT
    [x0,v0]=init_H_H(@modelequilibrium,dependentvar(ocEP),parametervalue(ocObj).',ap,varargin{:});
catch
    lasterr
    ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end