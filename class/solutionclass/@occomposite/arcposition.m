function out=arcposition(ocComp)

out=[];
if isempty(ocComp)
    return
end
out=cell(1,ocComp.order);
for ii=1:ocComp.order
    out{ii}=arcposition(ocComp.path{ii});
end