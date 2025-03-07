function ocObj=setparametervalue(ocObj,parval)

if isempty(ocObj)
    return
end
ocObj.Model=setparametervalue(ocObj.Model,parval);