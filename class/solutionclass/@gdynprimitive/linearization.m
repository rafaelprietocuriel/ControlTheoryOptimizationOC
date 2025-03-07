function out=linearization(gdynPrim)
%
out=linearization(gdynPrim.ocgtrajectory);

if isempty(out)
    ocObj=loadmodel(gdynPrim);
    if isempty(ocObj)
        out=[];
        return
    end
    par=modelparameter(gdynPrim);
    ocObj=changeparametervalue(ocObj,par);
    arcarg=arcargument(gdynPrim);
    depvar=dependentvar(gdynPrim);
    out{1}=feval(ocObj,'CanonicalSystemGeneralizedJacobian',0,depvar{1},par,arcarg(1));
end
