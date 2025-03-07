function arcarg=arcargument(dgObj)
%
% ARCARGUMENT

arcarg='';
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'argument');
arcarg=info.value;