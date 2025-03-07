function out=equalityconstraintnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation4staticoptmodel(ocObj.Model,'equalityconstraintnum');
out=info.value;