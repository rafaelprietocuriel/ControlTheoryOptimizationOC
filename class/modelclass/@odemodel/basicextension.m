function ext=basicextension(odeObj)

ext='';
if isempty(odeObj)
    return
end
ext=basicextension(odeObj.Model);
