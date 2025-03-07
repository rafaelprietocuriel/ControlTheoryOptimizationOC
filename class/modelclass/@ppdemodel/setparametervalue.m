function ppdeObj=setparametervalue(ppdeObj,parval)

if isempty(ppdeObj)
    return
end
ppdeObj.Model=setparametervalue(ppdeObj.Model,parval);