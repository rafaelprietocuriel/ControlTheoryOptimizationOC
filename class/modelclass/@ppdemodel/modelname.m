function name=modelname(ppdeocObj)

name='';
if isempty(ppdeocObj)
    return
end
name=modelname(ppdeocObj.Model);