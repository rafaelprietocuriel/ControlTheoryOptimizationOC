function [x0,v0]=ocmatinit_EPC_EPC(ocObj,dynPrim,ap,varargin)
%
% INIT_EP_EP initialization for equilibrium continuation
%
% [X0,V0]= INIT_EP_EP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT 


x0=[];
v0=[];
plotcoord=[];
arcarg=[];
if isempty(ocObj) || isempty(dynPrim)
    return
end
if ~isequilibrium(dynPrim)
    ocmaterror('Input argument is not an equilibrium.')
end
freeparameter=find(strcmpi(varargin,'freeparameter'),1);
plotcoordindex=find(strcmpi(varargin,'plotcoord'),1);
if ~isempty(freeparameter)
    freeparameterindex=parameterindex(ocObj,varargin{freeparameter+1});
end
if ~isempty(plotcoordindex)
    plotcoord=parameterindex(ocObj,varargin{plotcoordindex+1});
end

if ischar(ap)
    ap=parameterindex(ocObj,ap);
end
depvar=dependentvar(dynPrim);
parcoord=length(depvar);
if isempty(plotcoord)
    plotcoord=2;
end
freepar=parametervalue(ocObj,freeparameterindex);
depvar=[depvar;freepar];

OCMATFINITCONT.modeltype=modeltype(ocObj);
OCMATFINITCONT.parametervalue=parametervalue(ocObj);
OCMATFINITCONT.freeparameterindex=freeparameterindex;
OCMATFINITCONT.depvarcoord=1:parcoord;
OCMATFINITCONT.freeparametercoord=parcoord+1:length(depvar);
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
            [x0,v0]=init_EP_EP(@modelequilibriumc,depvar,parametervalue(ocObj),ap);
            %[x0,v0]=init_EP_EP(@modelequilibrium,depvar,parametervalue(ocObj),ap,varargin{3:end});
        case 'odemodel'
            [x0,v0]=init_EP_EP(@modelequilibriumc,depvar,parametervalue(ocObj),ap,varargin{3:end});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end