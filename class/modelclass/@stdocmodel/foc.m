function focStruct=foc(ocObj)

focStruct=[];
if isempty(ocObj)
    return
end
ocStruct=loadmodeldata(ocObj);
focStruct=ocStruct.foc;