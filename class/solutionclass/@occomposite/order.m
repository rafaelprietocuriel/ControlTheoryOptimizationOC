function out=order(ocComp)

out=[];
if isempty(ocComp)
    return
end
out=ocComp.order;