function out=localintegralstatenum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'localintegralstatenum');
out=info.value;