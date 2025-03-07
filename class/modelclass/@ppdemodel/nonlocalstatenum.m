function out=nonlocalstatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'nonlocalstatenum');
out=info.value;