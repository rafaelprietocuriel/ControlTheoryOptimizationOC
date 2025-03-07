function odeObjNew=copyresult(odeObjNew,odeObj)

if isempty(odeObj)
    odeObjNew=odeObj;
    return
end
if ~strcmp(modelname(odeObjNew),modelname(odeObj))
end
odeObjNew.Result=odeObj.Result;