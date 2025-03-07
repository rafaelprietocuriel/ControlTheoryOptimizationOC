function mtype=modeltype(ocObj)

mtype='';
if isempty(ocObj)
    return
end
mtype=modeltype(ocObj.Model);
