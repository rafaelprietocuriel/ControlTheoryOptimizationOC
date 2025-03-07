function varargout=cleanmodel(dgObj)

dgObj.Result=struct([]);

if ~nargout
    assignin('caller',inputname(1),dgObj);
else
    varargout{1}=dgObj;
end