function out=discountfactor(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'discountfactor');
out=info.value;