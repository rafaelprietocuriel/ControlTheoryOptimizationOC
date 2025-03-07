function out=constraintidentifier(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'inequalitycontrolconstraintidentifier');
out=info.value;