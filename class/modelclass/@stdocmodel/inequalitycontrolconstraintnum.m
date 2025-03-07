function out=inequalitycontrolconstraintnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'inequalitycontrolconstraintnum');
out=info.value;