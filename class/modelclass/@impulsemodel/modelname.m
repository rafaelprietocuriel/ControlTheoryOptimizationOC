function name=modelname(ocObj)

name='';
if isempty(ocObj)
    return
end
name=modelname(ocObj.Model);