function out=arcargument(ocComp)

out=[];
if isempty(ocComp)
    return
end
for ii=1:ocComp.order
    out{ii}=arcargument(ocComp.path{ii});
end