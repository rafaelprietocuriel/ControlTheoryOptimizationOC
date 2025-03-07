function [x0,v0]=initocmat_EPDE_EPDE(ocObj,ppdePrim,ap,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_EP_EP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];
idx=[];
userfunc='';
if isempty(ocObj) || isempty(ppdePrim)
    return
end
if ~isequilibrium(ppdePrim)
    ocmaterror('Input argument is not an equilibrium.')
end
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
if isempty(ap) || ap>parameternum(ocObj) || ap<1
    ocmaterror('Parameter index ''%s'' is not valid',char(ap))
end
if nargin>=4
    userfunc=varargin{1};
end
if nargin>=5
    idx=varargin{2};
end
depvar=dependentvar(ppdePrim);
if isempty(idx)
    idx=1:length(depvar);
end
depvar=depvar(idx);
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
if isempty(userfunc)
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquation');
            OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'EquilibriumEquationJacobian');
            OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'EquilibriumEquationParameterJacobian');
            OCMATFINITCONT.arcarg=arcargument(ppdePrim);
        case 'odemodel'
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'Dynamics');
            OCMATFINITCONT.jacobian=modelspecificfunc(ocObj,'DynamicsJacobian');
            OCMATFINITCONT.parameterjacobian=modelspecificfunc(ocObj,'DynamicsParameterJacobian');
            OCMATFINITCONT.arcarg=arcargument(ppdePrim);
    end
else
    funch=feval(modelspecificfunc(ocObj,userfunc));
    OCMATFINITCONT.dynamics=funch{1};
    OCMATFINITCONT.jacobian=funch{2};
    OCMATFINITCONT.parameterjacobian=funch{3};
    OCMATFINITCONT.arcarg=arcargument(ppdePrim);
end

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');
OCMATFINITCONT.plotcoord=2;
try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            [x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
        case 'odemodel'
            [x0,v0]=init_EP_EP(@odemodelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end