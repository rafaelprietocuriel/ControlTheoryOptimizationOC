function odeObj=setparametervalue(odeObj,parval)

if isempty(odeObj)
    return
end
odeObj.Model=setparametervalue(odeObj.Model,parval);