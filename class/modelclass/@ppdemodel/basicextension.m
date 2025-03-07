function ext=basicextension(ppdeObj)

ext='';
if isempty(ppdeObj)
    return
end
ext=basicextension(ppdeObj.Model);
