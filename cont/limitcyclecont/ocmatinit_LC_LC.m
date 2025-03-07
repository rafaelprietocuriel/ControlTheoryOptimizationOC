function [x0,v0]=ocmatinit_LC_LC(ocObj,dynPrim,ap,varargin)
%
% ocmatinit_LP_LP Initializes a fold continuation from a limitpoint.
%
% [X0,V0]= ocmatinit_LP_LP(OCOBJ,DYNPRIM,AP)
clear global OCMATFINITCONT 
global OCMATFINITCONT lds


x0=[];
v0=[];
if isempty(ocObj) || isempty(dynPrim)
    return
end
% if ~islimitcycle(dynPrim)
%     ocmaterror('Input argument is not an equilibrium.')
% end
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

modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');

try
    % calls MATCONT
    switch OCMATFINITCONT.modeltype
        case 'standardmodel'
            ntst=dynPrim.octrajectory.solverinfo.ntst;
            ncol=dynPrim.octrajectory.solverinfo.ncol;
            s.data.ntst=ntst;
            s.data.ncol=ncol;
            s.data.parametervalues=parametervalue(ocObj).';
            s.data.timemesh=dynPrim.octrajectory.solverinfo.tmesh;
            s.data.T=timehorizon(dynPrim);
            s.index=1;
            numcoeff=ntst*ncol*2*statenum(ocObj)+statenum(ocObj)*2+2;
            par=modelparameter(dynPrim);
            x=[dynPrim.octrajectory.y(:);period(dynPrim);par(ap)]; %dynPrim.octrajectory.solverinfo.coeff(1:numcoeff)
            [x0,v0]=init_LC_LC(@modellimitcycle,x,[],s,ap,ntst,ncol);
        case 'odemodel'
            [x0,v0]=init_LPC_LPC(@modellimitcycle,dependentvar(dynPrim),parametervalue(ocObj).',ap,varargin{:});
    end
catch
    lasterr
    %ocmaterror('MATCONT may not be installed or not on the MATLAB path.\n')
end