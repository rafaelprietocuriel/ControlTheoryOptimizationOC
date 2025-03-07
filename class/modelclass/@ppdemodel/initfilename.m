function name=initfilename(ppdeObj)

name='';
if isempty(ppdeObj)
    return
end
name=initfilename(ppdeObj.Model);