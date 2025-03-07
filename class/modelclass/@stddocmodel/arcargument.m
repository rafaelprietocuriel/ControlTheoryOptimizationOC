function arcarg=arcargument(ocObj)
%
% ARCARGUMENT

arcarg='';
if isempty(ocObj)
    return
end
info=retrievediffmodelinformation(ocObj.Model,'argument');
arcarg=info.value;