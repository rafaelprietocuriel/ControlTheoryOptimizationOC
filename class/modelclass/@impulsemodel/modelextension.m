function mext=modelextension(ocObj)

mext='';
if isempty(ocObj)
    return
end
mext=[basicextension(ocObj.Model) 'm'];
