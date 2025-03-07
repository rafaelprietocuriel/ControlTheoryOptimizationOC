function out=constraintnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation4staticoptmodel(ocObj.Model,'inequalityconstraintnum');
out=info.value;