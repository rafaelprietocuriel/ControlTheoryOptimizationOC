function mtype=modeltype(ocStruct)

mtype='';
if isempty(ocStruct)
    return
end
mtype=ocStruct.modeltype;