function varargout=cleanmodel(mmObj)

mmObj.Result=struct([]);

if ~nargout
    assignin('caller',inputname(1),mmObj);
else
    varargout{1}=mmObj;
end