function varargout=processinitfile(dgObj,varargin)

processinitfile(modelname(dgObj),varargin{:});
dgObj=stdocmodel(modelname(dgObj));
if ~nargout
    assignin('caller',inputname(1),dgObj);
else
    varargout{1}=dgObj;
end