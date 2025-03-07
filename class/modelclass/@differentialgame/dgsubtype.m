function subtype=dgsubtype(dgObj)

subtype='';
if isempty(dgObj)
    return
end
subtype=dgsubtype(dgObj.Model);