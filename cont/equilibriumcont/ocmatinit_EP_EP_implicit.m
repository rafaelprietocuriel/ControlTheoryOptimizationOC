function [x0,v0]=ocmatinit_EP_EP_implicit(dynPrim,ap,varargin)
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
plotcoord=[];
if isempty(dynPrim)
    return
end
if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
ocObj=loadmodel(dynPrim);
if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
OCMATFINITCONT.parametername=parametername(ocObj,ap);

if isempty(ap) || ap>parameternum(ocObj) || ap<1
    ocmaterror('Parameter index ''%s'' is not valid',char(ap))
end
if nargin>=3
    userfunc=varargin{1};
end
if nargin>=4
    idx=varargin{2};
end
if nargin>=5
    plotcoord=varargin{3};
end
if nargin>=6
    v0=varargin{4};
end
depvar=dependentvar(dynPrim);
depvar=depvar{1};
if isempty(idx)
    idx=1:length(depvar);
end
if isempty(plotcoord)
    plotcoord=1;
end
depvar=depvar(idx);
OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=modelparameter(dynPrim);
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
if isa(ap,'function_handle')
    OCMATFINITCONT.combinationfunction=ap;
    par=[parametervalue(ocObj) 0];
    ap=length(par);
    OCMATFINITCONT.parametername{1}='cont';

else
    OCMATFINITCONT.combinationfunction=[];
    par=parametervalue(ocObj);
    OCMATFINITCONT.parametername=parametername(ocObj,ap);

end
OCMATFINITCONT.activeparameter=ap;
try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case {'standardmodel','differentialgame'}
            x0=init_EP_EP(@modelequilibrium,depvar,modelparameter(dynPrim),ap);
            %[x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end