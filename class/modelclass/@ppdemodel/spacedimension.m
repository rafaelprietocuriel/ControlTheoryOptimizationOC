function spacedim=spacedimension(ppdeObj)
%
% SPACEDIMENSION

spacedim=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'spacedimension');
spacedim=info.value;