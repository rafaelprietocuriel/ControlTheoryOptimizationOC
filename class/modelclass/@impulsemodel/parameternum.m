function parnum=parameternum(ocObj)

parnum=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj,'parameternum');
parnum=info.value;