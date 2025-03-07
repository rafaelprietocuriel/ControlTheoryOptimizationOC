function ppdeObjNew=copyresult(ppdeObjNew,ppdeObj)

if isempty(ppdeObj)
    ppdeObjNew=ppdeObj;
    return
end
if ~strcmp(modelname(ppdeObjNew),modelname(ppdeObj))
end
ppdeObjNew.Result=ppdeObj.Result;