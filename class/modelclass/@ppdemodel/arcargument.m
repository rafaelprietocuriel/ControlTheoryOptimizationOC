function arcarg=arcargument(ppdeObj)
%
% ARCARGUMENT

arcarg='';
if isempty(ppdeObj)
    return
end
info=retrievemodelinformation(ppdeObj.Model,'argument');
arcarg=info.value;