function out=independentvar(ocComp)

out=[];
if isempty(ocComp)
    return
end
out=cell(1,ocComp.order);
for ii=1:ocComp.order
    out{ii}=independentvar(ocComp.path{ii});
end