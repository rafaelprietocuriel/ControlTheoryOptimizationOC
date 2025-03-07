function mtype=modelsubtype(dgObj)

mtype='';
if isempty(dgObj)
    return
end
subtype=retrievedifferentialgameinformation(dgObj.Model,'subtype');

mtype=subtype.value;
