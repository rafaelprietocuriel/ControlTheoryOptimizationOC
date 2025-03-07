function mext=modelextension(odeObj)

mext='';
if isempty(odeObj)
    return
end
mext=[basicextension(odeObj.Model) 'm'];
