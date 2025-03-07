function out=localcontrolnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'localcontrolnum');
out=info.value;