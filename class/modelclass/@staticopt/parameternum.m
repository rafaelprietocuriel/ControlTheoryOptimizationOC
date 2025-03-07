function parnum=parameternum(ocObj)

parnum=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'parameternum');
parnum=info.value;