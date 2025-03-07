function timearg=timeargument(ppdeObj)
%
% TIMEARGUMENT

timearg='';
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'time');
timearg=info.value;