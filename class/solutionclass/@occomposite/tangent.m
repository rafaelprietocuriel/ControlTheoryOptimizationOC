function out=tangent(ocComp)

out=[];
if isempty(ocComp)
    return
end
out=cell(1,ocComp.order);
for ii=1:ocComp.order
    out{ii}=tangent(ocComp.path{ii});
end