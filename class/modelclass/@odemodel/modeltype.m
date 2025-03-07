function mtype=modeltype(odeObj)

mtype='';
if isempty(odeObj)
    return
end
mtype=modeltype(odeObj.Model);
