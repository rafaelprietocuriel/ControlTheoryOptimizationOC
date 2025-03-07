function varargout=clearmodel(ocObj)

ocObj.Result=struct([]);

if ~nargout
    assignin('caller',inputname(1),ocObj);
else
    varargout{1}=ocObj;
end