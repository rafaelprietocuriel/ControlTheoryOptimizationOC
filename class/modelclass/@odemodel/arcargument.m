function arcarg=arcargument(odeObj)
%
% ARCARGUMENT

arcarg='';
if isempty(odeObj)
    return
end
info=retrieveodemodelinformation(odeObj.Model,'argument');
arcarg=info.value;