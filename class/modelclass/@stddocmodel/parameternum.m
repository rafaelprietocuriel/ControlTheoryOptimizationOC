function parnum=parameternum(ocObj)

parnum=[];
if isempty(ocObj)
    return
end
info=retrievediffmodelinformation(ocObj.Model,'parameternum');
parnum=info.value;