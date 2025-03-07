function varargout=feval(dynPrim,ocObj,basicfuncname,varargin)

if nargout==0
    feval(ocObj,basicfuncname,independentvar(dynPrim),dependentvar(dynPrim),parametervalue(ocObj),arcargument(dynPrim));
    try
        varargout{1}=ans;
    end
else
    [varargout{1:nargout}]=feval(ocObj,basicfuncname,independentvar(dynPrim),dependentvar(dynPrim),parametervalue(ocObj),arcargument(dynPrim));
end
