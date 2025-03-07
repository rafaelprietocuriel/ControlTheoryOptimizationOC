function out=localstatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'localstatenum');
out=info.value;