function out=nonlocalintegralstatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'nonlocalintegralstatenum');
out=info.value;