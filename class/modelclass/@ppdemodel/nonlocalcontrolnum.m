function out=nonlocalcontrolnum(ocObj)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'nonlocalcontrolnum');
out=info.value;