function ocObjNew=mergecontresult(ocObjNew,ocObj)

if isempty(ocObj)
    ocObjNew=ocObj;
    return
end
if ~strcmp(modelname(ocObjNew),modelname(ocObj))
end
ocObjNew.Result.Continuation=[ocObj.Result.Continuation ocObjNew.Result.Continuation];