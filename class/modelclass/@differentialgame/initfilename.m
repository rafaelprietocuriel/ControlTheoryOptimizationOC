function name=initfilename(dgObj)

name='';
if isempty(dgObj)
    return
end
name=initfilename(dgObj.Model);