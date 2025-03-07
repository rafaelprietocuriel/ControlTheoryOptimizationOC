function arcarg=arcargument(ocObj)
%
% ARCARGUMENT

arcarg='';
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'argument');
arcarg=info.value;