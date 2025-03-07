function out=controlnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrieveimpulsemodelinformation(ocObj.Model,'controlnum');
out=info.value;