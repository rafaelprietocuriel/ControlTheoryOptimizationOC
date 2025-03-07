function name=initfilename(ocObj)

name='';
if isempty(ocObj)
    return
end
name=initfilename(ocObj.Model);