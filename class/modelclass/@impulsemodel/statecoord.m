function out=statecoord(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj,'statenum');
out=1:info.value;