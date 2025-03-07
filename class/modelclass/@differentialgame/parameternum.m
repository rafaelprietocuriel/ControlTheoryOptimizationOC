function parnum=parameternum(dgObj)

parnum=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'parameternum');
parnum=info.value;