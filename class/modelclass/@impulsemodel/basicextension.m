function ext=basicextension(ocObj)

ext='';
if isempty(ocObj)
    return
end
ext=basicextension(ocObj.Model);
