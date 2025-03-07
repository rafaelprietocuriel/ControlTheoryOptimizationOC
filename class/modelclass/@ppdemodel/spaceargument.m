function spacearg=spaceargument(ppdeObj)
%
% SPACEARGUMENT

spacearg='';
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'space');
spacearg=info.value;