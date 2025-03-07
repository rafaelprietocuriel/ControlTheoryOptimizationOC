function ext=basicextension(dgObj)

ext='';
if isempty(dgObj)
    return
end
ext=basicextension(dgObj.Model);
