function ocObj=renameresult(ocObj,oldname,newname)
%
%
if ~isfield(ocObj.Result,oldname)
    ocmatmsg([oldname ' is not a name of the result field.'])
    return
end

if isfield(ocObj.Result,newname)
end
ocObj.Result.(newname)=ocObj.Result.(oldname);
ocObj.Result=rmfield(ocObj.Result,oldname);
