function out=costatecoord(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'statenum');
out=(1:info.value)+info.value;