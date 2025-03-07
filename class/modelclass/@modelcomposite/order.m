function out=order(mComp)

out=[];
if isempty(mComp)
    return
end
out=mComp.Order;