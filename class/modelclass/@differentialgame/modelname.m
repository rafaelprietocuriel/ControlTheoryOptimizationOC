function name=modelname(dgObj)

name='';
if isempty(dgObj)
    return
end
name=modelname(dgObj.Model);