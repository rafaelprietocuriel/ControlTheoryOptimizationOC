function out=impulsecontrolnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrieveimpulsemodelinformation(ocObj.Model,'icontrolnum');
out=info.value;