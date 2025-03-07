function parnum=parameternum(odeObj)

parnum=[];
if isempty(odeObj)
    return
end
info=retrieveodemodelinformation(odeObj.Model,'parameternum');
parnum=info.value;