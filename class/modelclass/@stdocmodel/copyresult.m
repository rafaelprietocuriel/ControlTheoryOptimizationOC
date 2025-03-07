function ocObjNew=copyresult(ocObjNew,ocObj)
% COPYRESULT copies the results of an object
%
% OCOBJNEW=COPYRESULT(OCOBJNEW,OCOBJ)
if isempty(ocObj)
    ocObjNew=ocObj;
    return
end
if ~strcmp(modelname(ocObjNew),modelname(ocObj))
end
ocObjNew.Result=ocObj.Result;