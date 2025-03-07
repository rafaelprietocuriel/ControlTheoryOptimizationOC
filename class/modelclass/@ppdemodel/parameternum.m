function parnum=parameternum(ppdeObj)

parnum=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'parameternum');
parnum=info.value;