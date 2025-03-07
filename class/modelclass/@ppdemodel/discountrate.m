function out=discountrate(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'discountrate');

out=parametervalue(ppdeObj,info.value);