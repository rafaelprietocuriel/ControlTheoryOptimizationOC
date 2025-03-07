function [x0,v0]=ocmatinit_EP_EMF(ocObj,dynPrim,acoord,varargin)
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
if isempty(ocObj) || isempty(dynPrim)
    return
end
if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
if isempty(acoord) || acoord>parameternum(ocObj) || acoord<1
    ocmaterror('Parameter index ''%s'' is not valid',char(acoord))
end
if nargin>=4
    userfunc=varargin{1};
end
if nargin>=5
    idx=varargin{2};
end
depvar=dependentvar(dynPrim);
if isempty(idx)
    idx=1:length(depvar);
end
depvar=depvar(idx);
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
OCMATFINITCONT.equilibriummfcoordposition=parameternum(ocObj)+1;
OCMATFINITCONT.equilibriummffreecoord=acoord;
OCMATFINITCONT.equilibriummfcoord=setdiff(idx,acoord);
OCMATFINITCONT.parametercoord=1:parameternum(ocObj);
totalpar=[parametervalue(ocObj) depvar(acoord)];
depvar(acoord)=[];
if isempty(userfunc)
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            OCMATFINITCONT.dynamics=modelspecificfunc(ocObj,'EquilibriumEquation');
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

try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            [x0,v0]=init_EP_EP(@modelequilibriummf,depvar,totalpar,parameternum(ocObj)+1,varargin{3:end});
        case 'odemodel'
            [x0,v0]=init_EP_EP(@odemodelequilibriummf,depvar,totalpar,acoord,varargin{3:end});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end