function out=nonlocalodestatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'nonlocalodestatenum');
out=info.value;