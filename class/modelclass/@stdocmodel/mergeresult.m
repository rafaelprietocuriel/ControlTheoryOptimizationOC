function ocObjNew=mergeresult(ocObjNew,ocObj)

if isempty(ocObj)
    ocObjNew=ocObj;
    return
end
if ~strcmp(modelname(ocObjNew),modelname(ocObj))
end
ocObjNew.Result=[ocObj.Result ocObjNew.Result];