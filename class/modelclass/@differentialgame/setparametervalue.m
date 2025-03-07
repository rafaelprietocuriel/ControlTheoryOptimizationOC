function dgObj=setparametervalue(dgObj,parval)

if isempty(dgObj)
    return
end
dgObj.Model=setparametervalue(dgObj.Model,parval);