function varargout=cleanmodel(odeObj)
odeObj.Result=struct([]);

if ~nargout
    assignin('caller',inputname(1),odeObj);
else
    varargout{1}=odeObj;
end