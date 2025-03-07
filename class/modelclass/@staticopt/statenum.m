function out=statenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation4staticoptmodel(ocObj.Model,'statenum');
out=info.value;