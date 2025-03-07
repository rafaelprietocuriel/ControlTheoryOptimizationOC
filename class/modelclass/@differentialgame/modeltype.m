function mtype=modeltype(dgObj)

mtype='';
if isempty(dgObj)
    return
end
mtype=modeltype(dgObj.Model);
