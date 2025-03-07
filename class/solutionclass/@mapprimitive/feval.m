function varargout=feval(mapPrim,ocObj,basicfuncname,varargin)

if nargout==0
    feval(ocObj,basicfuncname,independentvar(mapPrim),dependentvar(mapPrim),parametervalue(ocObj),arcargument(mapPrim));
    try
        varargout{1}=ans;
    end
else
    [varargout{1:nargout}]=feval(ocObj,basicfuncname,independentvar(mapPrim),dependentvar(mapPrim),parametervalue(ocObj),arcargument(mapPrim));
    end
