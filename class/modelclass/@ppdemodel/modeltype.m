function mtype=modeltype(ppdeocObj)

mtype='';
if isempty(ppdeocObj)
    return
end
mtype=modeltype(ppdeocObj.Model);