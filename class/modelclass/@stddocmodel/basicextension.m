function ext=basicextension(docObj)

ext='';
if isempty(docObj)
    return
end
ext=basicextension(docObj.Model);
