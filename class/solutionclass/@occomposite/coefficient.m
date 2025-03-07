function out=coefficient(ocComp)

out=[];
if isempty(ocComp)
    return
end
out=cell(1,ocComp.order);
for ii=1:ocComp.order
    out{ii}=coefficient(ocComp.path{ii});
end