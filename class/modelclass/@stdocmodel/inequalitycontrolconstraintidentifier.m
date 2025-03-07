function out=inequalitycontrolconstraintidentifier(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'inequalitycontrolconstraintidentifier');
out=info.value;