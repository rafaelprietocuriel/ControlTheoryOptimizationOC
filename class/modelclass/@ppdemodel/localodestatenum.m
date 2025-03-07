function out=localodestatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'localodestatenum');
out=info.value;