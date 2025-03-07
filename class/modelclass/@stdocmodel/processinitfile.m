function varargout=processinitfile(ocObj,varargin)

processinitfile(modelname(ocObj),varargin{:});
ocObj=stdocmodel(modelname(ocObj));
if ~nargout
    assignin('caller',inputname(1),ocObj);
else
    varargout{1}=ocObj;
end