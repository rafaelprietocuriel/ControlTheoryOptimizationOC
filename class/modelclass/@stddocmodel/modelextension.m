function mext=modelextension(docObj)

mext='';
if isempty(docObj)
    return
end
mext=[basicextension(docObj.Model) 'm'];
